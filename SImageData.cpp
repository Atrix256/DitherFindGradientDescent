#define _CRT_SECURE_NO_WARNINGS

#include "SImageData.h"

#include <windows.h>  // for bitmap headers.  Sorry non windows people!
#include <thread>
#include <atomic>
#include <array>
#include "Misc.h"
#include "GaussianBlurGeneric.h"

const float c_pi = 3.14159265359f;

//======================================================================================
std::complex<float> DFTPixel (const SImageData &srcImage, size_t K, size_t L)
{
    std::complex<float> ret(0.0f, 0.0f);
  
    for (size_t x = 0; x < srcImage.m_width; ++x)
    {
        for (size_t y = 0; y < srcImage.m_height; ++y)
        {
            // Get the pixel value (assuming greyscale) and convert it to [0,1] space
            const uint8 *src = &srcImage.m_pixels[(y * srcImage.m_pitch) + x * 3];
            float grey = float(src[0]) / 255.0f;
  
            // Add to the sum of the return value
            float v = float(K * x) / float(srcImage.m_width);
            v += float(L * y) / float(srcImage.m_height);
            ret += std::complex<float>(grey, 0.0f) * std::polar<float>(1.0f, -2.0f * c_pi * v);
        }
    }
  
    return ret;
}
  
//======================================================================================
void ImageDFT (const SImageData &srcImage, SImageDataComplex &destImage)
{
    // NOTE: this function assumes srcImage is greyscale, so works on only the red component of srcImage.
    // ImageToGrey() will convert an image to greyscale.
 
    // size the output dft data
    destImage.m_width = srcImage.m_width;
    destImage.m_height = srcImage.m_height;
    destImage.m_pixels.resize(destImage.m_width*destImage.m_height);
 
    size_t numThreads = std::thread::hardware_concurrency();
    //if (numThreads > 0)
        //numThreads = numThreads - 1;
 
    std::vector<std::thread> threads;
    threads.resize(numThreads);
 
    printf("Doing DFT with %zu threads...\n", numThreads);
 
    // calculate 2d dft (brute force, not using fast fourier transform) multithreadedly
    std::atomic<size_t> nextRow(0);
    for (std::thread& t : threads)
    {
        t = std::thread(
            [&] ()
            {
                size_t row = nextRow.fetch_add(1);
                bool reportProgress = (row == 0);
                int lastPercent = -1;
 
                while (row < srcImage.m_height)
                {
                    // calculate the DFT for every pixel / frequency in this row
                    for (size_t x = 0; x < srcImage.m_width; ++x)
                    {
                        destImage.m_pixels[row * destImage.m_width + x] = DFTPixel(srcImage, x, row);
                    }
 
                    // report progress if we should
                    if (reportProgress)
                    {
                        int percent = int(100.0f * float(row) / float(srcImage.m_height));
                        if (lastPercent != percent)
                        {
                            lastPercent = percent;
                            printf("            \rDFT: %i%%", lastPercent);
                        }
                    }
 
                    // go to the next row
                    row = nextRow.fetch_add(1);
                }
            }
        );
    }
 
    for (std::thread& t : threads)
        t.join();
 
    printf("\n");
}
 
//======================================================================================
void GetMagnitudeData (const SImageDataComplex& srcImage, SImageData& destImage)
{
    // size the output image
    destImage.m_width = srcImage.m_width;
    destImage.m_height = srcImage.m_height;
    destImage.m_pitch = 4 * ((srcImage.m_width * 24 + 31) / 32);
    destImage.m_pixels.resize(destImage.m_pitch*destImage.m_height);
  
    // get floating point magnitude data
    std::vector<float> magArray;
    magArray.resize(srcImage.m_width*srcImage.m_height);
    float maxmag = 0.0f;
    for (size_t x = 0; x < srcImage.m_width; ++x)
    {
        for (size_t y = 0; y < srcImage.m_height; ++y)
        {
            // Offset the information by half width & height in the positive direction.
            // This makes frequency 0 (DC) be at the image origin, like most diagrams show it.
            int k = (x + (int)srcImage.m_width / 2) % (int)srcImage.m_width;
            int l = (y + (int)srcImage.m_height / 2) % (int)srcImage.m_height;
            const std::complex<float> &src = srcImage.m_pixels[l*srcImage.m_width + k];
  
            float mag = std::abs(src);
            if (mag > maxmag)
                maxmag = mag;
  
            magArray[y*srcImage.m_width + x] = mag;
        }
    }
    if (maxmag == 0.0f)
        maxmag = 1.0f;
  
    const float c = 255.0f / log(1.0f+maxmag);
  
    // normalize the magnitude data and send it back in [0, 255]
    for (size_t x = 0; x < srcImage.m_width; ++x)
    {
        for (size_t y = 0; y < srcImage.m_height; ++y)
        {
            float src = c * log(1.0f + magArray[y*srcImage.m_width + x]);
  
            uint8 magu8 = uint8(src);
  
            uint8* dest = &destImage.m_pixels[y*destImage.m_pitch + x * 3];
            dest[0] = magu8;
            dest[1] = magu8;
            dest[2] = magu8;
        }
    }
}

//======================================================================================
bool ImageSave (const SImageData &image, const char *fileName)
{
    // open the file if we can
    FILE *file;
    file = fopen(fileName, "wb");
    if (!file) {
        printf("Could not save %s\n", fileName);
        return false;
    }
   
    // make the header info
    BITMAPFILEHEADER header;
    BITMAPINFOHEADER infoHeader;
   
    header.bfType = 0x4D42;
    header.bfReserved1 = 0;
    header.bfReserved2 = 0;
    header.bfOffBits = 54;
   
    infoHeader.biSize = 40;
    infoHeader.biWidth = (LONG)image.m_width;
    infoHeader.biHeight = (LONG)image.m_height;
    infoHeader.biPlanes = 1;
    infoHeader.biBitCount = 24;
    infoHeader.biCompression = 0;
    infoHeader.biSizeImage = (DWORD) image.m_pixels.size();
    infoHeader.biXPelsPerMeter = 0;
    infoHeader.biYPelsPerMeter = 0;
    infoHeader.biClrUsed = 0;
    infoHeader.biClrImportant = 0;
   
    header.bfSize = infoHeader.biSizeImage + header.bfOffBits;
   
    // write the data and close the file
    fwrite(&header, sizeof(header), 1, file);
    fwrite(&infoHeader, sizeof(infoHeader), 1, file);
    fwrite(&image.m_pixels[0], infoHeader.biSizeImage, 1, file);
    fclose(file);
  
    return true;
}

//======================================================================================
bool ImageLoad (const char *fileName, SImageData& imageData)
{
    // open the file if we can
    FILE *file;
    file = fopen(fileName, "rb");
    if (!file)
        return false;
 
    // read the headers if we can
    BITMAPFILEHEADER header;
    BITMAPINFOHEADER infoHeader;
    if (fread(&header, sizeof(header), 1, file) != 1 ||
        fread(&infoHeader, sizeof(infoHeader), 1, file) != 1 ||
        header.bfType != 0x4D42 || infoHeader.biBitCount != 24)
    {
        fclose(file);
        return false;
    }
 
    // read in our pixel data if we can. Note that it's in BGR order, and width is padded to the next power of 4
    imageData.m_pixels.resize(infoHeader.biSizeImage);
    fseek(file, header.bfOffBits, SEEK_SET);
    if (fread(&imageData.m_pixels[0], imageData.m_pixels.size(), 1, file) != 1)
    {
        fclose(file);
        return false;
    }
 
    imageData.m_width = infoHeader.biWidth;
    imageData.m_height = infoHeader.biHeight;
    imageData.m_pitch = 4 * ((imageData.m_width * 24 + 31) / 32);
 
    fclose(file);
    return true;
}

//======================================================================================
void ImageInit (SImageData& image, size_t width, size_t height)
{
    image.m_width = width;
    image.m_height = height;
    image.m_pitch = 4 * ((width * 24 + 31) / 32);
    image.m_pixels.resize(image.m_pitch * image.m_height);
    std::fill(image.m_pixels.begin(), image.m_pixels.end(), 0);
}

//======================================================================================
void ImageConvertToLuma (SImageData& image)
{
    ImageForEachPixel(
        image,
        [] (SColor& pixel, size_t pixelIndex, size_t x, size_t y)
        {
            uint8 lumau8 = pixel.LumaU8();
            pixel.R = lumau8;
            pixel.G = lumau8;
            pixel.B = lumau8;
        }
    );
}

//======================================================================================
void ImageCombine2 (const SImageData& imageA, const SImageData& imageB, SImageData& result)
{
    // put the images side by side. A on left, B on right
    ImageInit(result, imageA.m_width + imageB.m_width, max(imageA.m_height, imageB.m_height));
    std::fill(result.m_pixels.begin(), result.m_pixels.end(), 0);

    // image A on left
    for (size_t y = 0; y < imageA.m_height; ++y)
    {
        SColor* destPixel = (SColor*)&result.m_pixels[y * result.m_pitch];
        SColor* srcPixel = (SColor*)&imageA.m_pixels[y * imageA.m_pitch];
        for (size_t x = 0; x < imageA.m_width; ++x)
        {
            destPixel[0] = srcPixel[0];
            ++destPixel;
            ++srcPixel;
        }
    }

    // image B on right
    for (size_t y = 0; y < imageB.m_height; ++y)
    {
        SColor* destPixel = (SColor*)&result.m_pixels[y * result.m_pitch + imageA.m_width * 3];
        SColor* srcPixel = (SColor*)&imageB.m_pixels[y * imageB.m_pitch];
        for (size_t x = 0; x < imageB.m_width; ++x)
        {
            destPixel[0] = srcPixel[0];
            ++destPixel;
            ++srcPixel;
        }
    }
}

//======================================================================================
void ImageCombine3 (const SImageData& imageA, const SImageData& imageB, const SImageData& imageC, SImageData& result)
{
    SImageData temp;
    ImageCombine2(imageA, imageB, temp);
    ImageCombine2(temp, imageC, result);
}

//======================================================================================
void ImageCombine5 (const SImageData& imageA, const SImageData& imageB, const SImageData& imageC, const SImageData& imageD, const SImageData& imageE, SImageData& result)
{
    SImageData temp;
    ImageCombine2(imageA, imageB, temp);
    ImageCombine2(temp, imageC, result);
    ImageCombine2(result, imageD, temp);
    ImageCombine2(temp, imageE, result);
}

//======================================================================================
void ImageDither (const SImageData& sourceImage, const SImageData& noiseImage, SImageData& result)
{
    // init the result image
    ImageInit(result, sourceImage.m_width, sourceImage.m_height);

    // make the result image
    for (size_t y = 0; y < sourceImage.m_height; ++y)
    {
        SColor* srcDitherPixel = (SColor*)&sourceImage.m_pixels[y * sourceImage.m_pitch];
        SColor* destDitherPixel = (SColor*)&result.m_pixels[y * result.m_pitch];

        for (size_t x = 0; x < sourceImage.m_width; ++x)
        {
            // tile the noise in case it isn't the same size as the image we are dithering
            size_t noiseX = x % noiseImage.m_width;
            size_t noiseY = y % noiseImage.m_height;
            SColor* noisePixel = (SColor*)&noiseImage.m_pixels[noiseY * noiseImage.m_pitch + noiseX * 3];

            uint8 value = 0;
            if (noisePixel->R < srcDitherPixel->R)
                value = 255;

            destDitherPixel->R = value;
            destDitherPixel->G = value;
            destDitherPixel->B = value;

            ++srcDitherPixel;
            ++destDitherPixel;
        }
    }
}

//======================================================================================
void ImageGaussianBlur (const SImageData& sourceImage, SImageData& destImage, float blurSigma)
{
    // initialize a temporary and destination image
    SImageData tempImage;
    ImageInit(tempImage, sourceImage.m_width, sourceImage.m_height);
    ImageInit(destImage, sourceImage.m_width, sourceImage.m_height);

    // blur!
    GaussianBlurGeneric<SImageData, std::array<float, 3>>(sourceImage.m_width, sourceImage.m_height, blurSigma, sourceImage, tempImage, destImage,
        
        // AccumulateBlur
        [] (std::array<float, 3>& blurredPixel, const SImageData& sourceImage, int pixelX, int pixelY, float weight, bool firstCall)
        {
            if (firstCall)
                blurredPixel[0] = blurredPixel[1] = blurredPixel[2] = 0.0f;

            const SColor& pixel = sourceImage.GetPixel(pixelX, pixelY);
            blurredPixel[0] += float(pixel.R) * weight;
            blurredPixel[1] += float(pixel.G) * weight;
            blurredPixel[2] += float(pixel.B) * weight;
        },

        // WriteBlurredPixel
        [] (const std::array<float, 3>& blurredPixel, SImageData& destImage, int pixelX, int pixelY)
        {
            SColor& destPixel = destImage.GetPixel(pixelX, pixelY);
            destPixel.R = uint8(blurredPixel[0]);
            destPixel.G = uint8(blurredPixel[1]);
            destPixel.B = uint8(blurredPixel[2]);
        }
    );
}