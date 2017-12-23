#define _CRT_SECURE_NO_WARNINGS

#include "SImageData.h"
#include "CDualNumber.h"
#include "Misc.h"
#include "GaussianBlurGeneric.h"
#include <random>
#include <atomic>

const float c_pi = 3.14159265359f;

template <size_t N>
void MakeDitherImage (SImageData& ditherImage, const std::array<CDualNumber<N*N>, N*N>& ditherPattern, size_t width, size_t height)
{
    ImageInit(ditherImage, width, height);

    ImageForEachPixel(
        ditherImage,
        [&] (SColor& pixel, size_t pixelIndex, size_t x, size_t y)
        {
            size_t ditherX = x % N;
            size_t ditherY = y % N;
            uint8 ditherValue = uint8(ditherPattern[ditherY*N + ditherX].m_value * 255.0f);
            pixel.Set(ditherValue, ditherValue, ditherValue);
        }
    );
}

float CalculateActualMSE (const SImageData& sourceImageBlurred, const SImageData& ditheredImageBlurred)
{
    const size_t c_numPixels = sourceImageBlurred.m_width * sourceImageBlurred.m_height;

    float meanSquaredError = 0.0f;
    for (size_t pixelIndex = 0; pixelIndex < c_numPixels; ++pixelIndex)
    {
        // get the source pixel's luma from 0 to 1. The texture is already greyscale so luma is R
        float error = float(ditheredImageBlurred.GetPixel(pixelIndex).R)/255.0f - float(sourceImageBlurred.GetPixel(pixelIndex).R)/255.0f;
        float errorSquared = error * error;
        meanSquaredError = Lerp(meanSquaredError, errorSquared, 1.0f / float(pixelIndex + 1));
    }
    return meanSquaredError;
}

template <size_t N>
CDualNumber<N*N> CalculateMSE (const SImageData& sourceImage, const SImageData& sourceImageBlurred, const std::array<CDualNumber<N*N>, N*N>& ditherPattern, float blurSigma, float percent)
{
    typedef CDualNumber<N*N> TDualNumber;
    const size_t c_numPixels = sourceImage.m_width * sourceImage.m_height;

    // dither the source image
    std::vector<TDualNumber> ditheredResults;
    ditheredResults.resize(c_numPixels);
    for (size_t i = 0; i < c_numPixels; ++i)
    {
        // get the source pixel's luma from 0 to 1. The texture is already greyscale so luma is R
        size_t pixelX = i % sourceImage.m_height;
        size_t pixelY = i / sourceImage.m_height;
        TDualNumber source(float(sourceImage.GetPixel(pixelX, pixelY).R) / 255.0f);

        // get the dither pattern value
        size_t ditherX = pixelX % N;
        size_t ditherY = pixelY % N;
        TDualNumber dither = ditherPattern[ditherY * N + ditherX];

        // Doing a differentiable version of this dithering logic:
        // result = dither < src ? 1.0 : 0.0
        // 0.5 + atan(100*x) / pi  ====> pretty close to the step function, but still differentiable.
        // If you make the multiplier of x larger, it gets closer to the step function, but the derivatives get flatter which seems not so great since it would slow down gradient descent.
        // Note: we are starting the multiplier at 1 and going up to 100 so that it starts out loose and gets tighter. This is akin to "Graduated Optimization"
        // https://en.m.wikipedia.org/wiki/Graduated_optimization
        //TDualNumber multiplier = Lerp(TDualNumber(1.0f), TDualNumber(10000.0f), TDualNumber(percent));
        TDualNumber multiplier = TDualNumber(10000.0f);
        ditheredResults[i] = TDualNumber(0.5f) + atan(multiplier * (source - dither)) / TDualNumber(c_pi);
    }

    // blur the dithered image
    std::vector<TDualNumber> ditheredResultsBlurredTemp;
    std::vector<TDualNumber> ditheredResultsBlurred;
    ditheredResultsBlurredTemp.resize(c_numPixels);
    ditheredResultsBlurred.resize(c_numPixels);
    GaussianBlurGeneric<std::vector<TDualNumber>, TDualNumber>(sourceImage.m_width, sourceImage.m_height, blurSigma, ditheredResults, ditheredResultsBlurredTemp, ditheredResultsBlurred,

        // AccumulateBlur
        [&sourceImage] (TDualNumber& blurredPixel, const std::vector<TDualNumber>& source, int pixelX, int pixelY, float weight, bool firstCall)
        {
            if (firstCall)
                blurredPixel = TDualNumber(0.0f);

            const TDualNumber& sourcePixel = source[pixelY * sourceImage.m_width + pixelX];
            blurredPixel = blurredPixel + sourcePixel * TDualNumber(weight);
        },

        // WriteBlurredPixel
        [&sourceImage] (const TDualNumber& blurredPixel, std::vector<TDualNumber>& dest, int pixelX, int pixelY)
        {
            TDualNumber& destPixel = dest[pixelY * sourceImage.m_width + pixelX];
            destPixel = blurredPixel;
        }
    );

    // calculate the MSE of the blurred dither image, vs the blurred source image
    TDualNumber MSE = 0.0f;
    for (size_t i = 0; i < c_numPixels; ++i)
    {
        const TDualNumber& result = ditheredResultsBlurred[i];
        TDualNumber source = TDualNumber(float(sourceImageBlurred.GetPixel(i).R) / 255.0f);

        // calculate error squared
        TDualNumber error = result - source;
        TDualNumber errorSquared = error * error;

        // incrementally average MSE
        MSE = Lerp(MSE, errorSquared, TDualNumber(1.0f / float(i + 1)));
    }

    return MSE;
}

template <size_t N>
void SaveState (const SImageData& sourceImage, const SImageData& sourceImageBlurred, const std::array<CDualNumber<N*N>, N*N>& ditherPattern, float blurSigma, const char* imageFileName, const char* textFileName)
{
    typedef CDualNumber<N*N> TDualNumber;

    // make the dither pattern into an image
    SImageData ditherImage;
    MakeDitherImage<N>(ditherImage, ditherPattern, sourceImage.m_width, sourceImage.m_height);

    // use the dither pattern image to dither the source
    SImageData sourceImageDithered;
    ImageDither(sourceImage, ditherImage, sourceImageDithered);

    // blur the dithered image
    SImageData sourceImageDitheredBlurred;
    ImageGaussianBlur(sourceImageDithered, sourceImageDitheredBlurred, blurSigma);

    // make the output image and save it
    SImageData outputImage;
    ImageCombine5(sourceImage, ditherImage, sourceImageDithered, sourceImageDitheredBlurred, sourceImageBlurred, outputImage);
    ImageSave(outputImage, imageFileName);

    // calculate the MSE and actual MSE
    TDualNumber meanSquaredError = CalculateMSE<N>(sourceImage, sourceImageBlurred, ditherPattern, blurSigma, 1.0f);
    float actualMSE = CalculateActualMSE(sourceImageBlurred, sourceImageDitheredBlurred);

    // write MSE values and pattern
    FILE *file = fopen(textFileName, "w+t");
    fprintf(file, "Branchless (Atan Step) MSE = %f\n", meanSquaredError.m_value);
    fprintf(file, "Actual MSE = %f\n", actualMSE);
    for (size_t y = 0; y < N; ++y)
    {
        for (size_t x = 0; x < N; ++x)
        {
            fprintf(file, "%f ", ditherPattern[y * N + x].m_value);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

template <size_t N>
void FindDitherPattern (const SImageData& sourceImage, size_t numIterations, float learningRateStart, float learningRateEnd, float blurSigma)
{
    typedef CDualNumber<N*N> TDualNumber;
    const size_t c_numPixels = sourceImage.m_width * sourceImage.m_height;

    // blur the source image
    SImageData sourceImageBlurred;
    ImageInit(sourceImageBlurred, sourceImage.m_width, sourceImage.m_height);
    ImageGaussianBlur(sourceImage, sourceImageBlurred, blurSigma);

    // initialize the dither pattern to random numbers, but set the derivative to 1 for whichever pixel it is
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    std::array<TDualNumber, N*N> ditherPattern;
    for (size_t i = 0; i < N*N; ++i)
    {
        TDualNumber& dualNumber = ditherPattern[i];
        dualNumber.m_value = dist(rng);
        dualNumber.m_derivatives[i] = 1.0f;
    }

    // calculate and report starting state
    printf("Starting Branchless (Atan Step) MSE = %f\n", CalculateMSE<N>(sourceImage, sourceImageBlurred, ditherPattern, blurSigma, 1.0f).m_value);

    // save our starting point
    SaveState<N>(sourceImage, sourceImageBlurred, ditherPattern, blurSigma, "start.bmp", "start.txt");

    // calculate the starting MSE
    TDualNumber meanSquaredError = CalculateMSE<N>(sourceImage, sourceImageBlurred, ditherPattern, blurSigma, 0.0f);

    // store off the current pattern as the best pattern
    std::array<TDualNumber, N*N> bestDitherPattern = ditherPattern;
    float bestDitherPatternMSE = meanSquaredError.m_value;
    size_t bestDitherPatternIndex = 0;

    // do the specified number of iterations
    for (size_t i = 0; i < numIterations; ++i)
    {
        // calculate the learning learning rate
        float percent = float(i + 1) / float(numIterations);
        float learningRate = Lerp(learningRateStart, learningRateEnd, percent);

        // Apply gradient descent to the dither pattern
        for (size_t derivativeIndex = 0; derivativeIndex < N*N; ++derivativeIndex)
        {
            ditherPattern[derivativeIndex].m_value -= meanSquaredError.m_derivatives[derivativeIndex] * learningRate;
            ditherPattern[derivativeIndex].m_value = Clamp(ditherPattern[derivativeIndex].m_value, 0.0f, 1.0f);
        }

        // Calculate and report current state
        meanSquaredError = CalculateMSE<N>(sourceImage, sourceImageBlurred, ditherPattern, blurSigma, percent);
        printf("\r[%i%%] Branchless (Atan Step) MSE = %f", int(100.0f * percent), meanSquaredError.m_value);

        // keep track of the best pattern we've seen
        if (meanSquaredError.m_value < bestDitherPatternMSE)
        {
            bestDitherPatternMSE = meanSquaredError.m_value;
            bestDitherPatternIndex = i;
            bestDitherPattern = ditherPattern;
        }
    }

    // report the dither pattern we are using
    printf("\n\nBest dither pattern found at percent %i%%.\nBranchless (Atan Step) MSE = %f\n\n", int(100.0f *float(bestDitherPatternIndex + 1) / float(numIterations)), bestDitherPatternMSE);
    ditherPattern = bestDitherPattern;

    // save our result
    SaveState<N>(sourceImage, sourceImageBlurred, ditherPattern, blurSigma, "stop.bmp", "stop.txt");
}

int main (int argc, char** argv)
{
    SImageData sourceImage;
    if (!ImageLoad("src/ditherimage.bmp", sourceImage))
    {
        printf("Could not load src/ditherimage.bmp");
        return 0;
    }
    ImageConvertToLuma(sourceImage);

    FindDitherPattern<3>(sourceImage, 1000, 3.0f, 0.1f, 2.0f);

    system("Pause");

    return 0;
}