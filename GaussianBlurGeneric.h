// A generic Gaussian blur.  Needed so it can be used to blur SImageData, but also so it can blur dual number pixels and be differentiable!
#pragma once

#include "Misc.h"

#include <cmath>
#include <vector>
#include <atomic>

namespace GaussianBlurGenericPrivate
{
    inline float Gaussian(float sigma, float x)
    {
        return expf(-(x*x) / (2.0f * sigma*sigma));
    }

    inline float GaussianSimpsonIntegration(float sigma, float a, float b)
    {
        return
            ((b - a) / 6.0f) *
            (Gaussian(sigma, a) + 4.0f * Gaussian(sigma, (a + b) / 2.0f) + Gaussian(sigma, b));
    }

    inline std::vector<float> GaussianKernelIntegrals(float sigma, int taps)
    {
        std::vector<float> ret;
        float total = 0.0f;
        for (int i = 0; i < taps; ++i)
        {
            float x = float(i) - float(taps / 2);
            float value = GaussianSimpsonIntegration(sigma, x - 0.5f, x + 0.5f);
            ret.push_back(value);
            total += value;
        }
        // normalize it
        for (unsigned int i = 0; i < ret.size(); ++i)
        {
            ret[i] /= total;
        }
        return ret;
    }

    template <typename TPIXELSTORAGE, typename TBLURACCUMULATE, typename TACCUMULATEBLURLAMBDA, typename TWRITEBLURREDPIXELLAMBDA, bool X_AXIS>
    static void BlurOneAxis (const TPIXELSTORAGE& sourceImage, TPIXELSTORAGE& destImage, size_t width, size_t height, float blurSigma, int blurSize, TACCUMULATEBLURLAMBDA& AccumulateBlur, TWRITEBLURREDPIXELLAMBDA& WriteBlurredPixel)
    {
        // calculate gaussian values
        std::vector<float> row = GaussianKernelIntegrals(blurSigma, blurSize);
        const int c_startOffset = -1 * int(row.size() / 2);

        // calculate the result for each pixel in the destination image, dividing the rows up across threads
        std::atomic<size_t> nextY(0);
        ForkJoin(
            [&nextY, &row, &AccumulateBlur, &WriteBlurredPixel, &sourceImage, &destImage, width, height, c_startOffset] ()
            {
                int y = (int)nextY.fetch_add(1);
                while (y < (int)height)
                {
                    for (int x = 0; x < (int)width; ++x)
                    {
                        // blur the source pixels to get the destination pixel value
                        TBLURACCUMULATE blurredPixel;
                        for (unsigned int i = 0; i < row.size(); ++i)
                        {
                            int pixelX = (X_AXIS == true) ? Clamp<int>(x + c_startOffset + i, 0, (int)width - 1) : x;
                            int pixelY = (X_AXIS != true) ? Clamp<int>(y + c_startOffset + i, 0, (int)height - 1) : y;

                            AccumulateBlur(blurredPixel, sourceImage, pixelX, pixelY, row[i], i == 0);
                        }

                        // write the destination pixel
                        WriteBlurredPixel(blurredPixel, destImage, x, y);
                    }
                    y = (int)nextY.fetch_add(1);
                }
            }
        );
    }
};

//======================================================================================
template <typename TPIXELSTORAGE, typename TBLURACCUMULATE, typename TACCUMULATEBLURLAMBDA, typename TWRITEBLURREDPIXELLAMBDA>
void GaussianBlurGeneric (size_t width, size_t height, float blurSigma, const TPIXELSTORAGE& source, TPIXELSTORAGE& temp, TPIXELSTORAGE& dest, TACCUMULATEBLURLAMBDA& AccumulateBlur, TWRITEBLURREDPIXELLAMBDA& WriteBlurredPixel)
{
    // no blur if zero or negative
    if (blurSigma <= 0.0f)
    {
        dest = source;
        return;
    }

    // calculate the number of pixels needed to represent a gaussian kernal that has values
    // down to the threshold amount.  A gaussian function technically has values everywhere
    // on the image, but the threshold lets us cut it off where the pixels contribute to
    // only small amounts that aren't as noticeable.
    // Also make sure the blur size is odd.
    const float c_threshold = 0.005f; // 0.5%
    int blurSize = int(std::floor(1.0f + 2.0f * std::sqrtf(-2.0f * blurSigma * blurSigma * std::log(c_threshold)))) + 1;
    blurSize = blurSize | 1;

    // horizontal blur from source into temp
    GaussianBlurGenericPrivate::BlurOneAxis<TPIXELSTORAGE, TBLURACCUMULATE, TACCUMULATEBLURLAMBDA, TWRITEBLURREDPIXELLAMBDA, true>(source, temp, width, height, blurSigma, blurSize, AccumulateBlur, WriteBlurredPixel);

    // vertical blur from temp into dest
    GaussianBlurGenericPrivate::BlurOneAxis<TPIXELSTORAGE, TBLURACCUMULATE, TACCUMULATEBLURLAMBDA, TWRITEBLURREDPIXELLAMBDA, false>(temp, dest, width, height, blurSigma, blurSize, AccumulateBlur, WriteBlurredPixel);
}