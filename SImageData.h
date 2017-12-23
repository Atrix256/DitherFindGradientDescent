#pragma once

#include <vector>
#include <stdint.h>
#include <complex>

//======================================================================================

typedef uint8_t uint8;

struct SColor;

//======================================================================================
struct SImageData
{
    SImageData ()
        : m_width(0)
        , m_height(0)
    { }

    SColor& GetPixel (size_t x, size_t y)
    {
        return *(SColor*)&m_pixels[y * m_pitch + x * 3];
    }

    const SColor& GetPixel (size_t x, size_t y) const
    {
        return *(SColor*)&m_pixels[y * m_pitch + x * 3];
    }

    SColor& GetPixel (size_t index)
    {
        size_t y = index / m_width;
        size_t x = index % m_width;
        return GetPixel(x, y);
    }

    const SColor& GetPixel (size_t index) const
    {
        size_t y = index / m_width;
        size_t x = index % m_width;
        return GetPixel(x, y);
    }
   
    size_t m_width;
    size_t m_height;
    size_t m_pitch;
    std::vector<uint8> m_pixels;
};
 
//======================================================================================
struct SColor
{
    SColor (uint8 _R = 0, uint8 _G = 0, uint8 _B = 0)
        : R(_R), G(_G), B(_B)
    { }

    inline void Set (uint8 _R, uint8 _G, uint8 _B)
    {
        R = _R;
        G = _G;
        B = _B;
    }

    float LumaFloat() const
    {
        return float(R) * 0.3f + float(G) * 0.59f + float(B) * 0.11f;
    }

    uint8 LumaU8 () const
    {
        return uint8(LumaFloat() + 0.5f);
    }
 
    uint8 B, G, R;
};

//======================================================================================
struct SImageDataComplex
{
    SImageDataComplex ()
        : m_width(0)
        , m_height(0)
    { }
  
    size_t m_width;
    size_t m_height;
    std::vector<std::complex<float>> m_pixels;
};

//======================================================================================
template <typename LAMBDA>
void ImageForEachPixel (SImageData& image, const LAMBDA& lambda)
{
    size_t pixelIndex = 0;
    for (size_t y = 0; y < image.m_height; ++y)
    {
        SColor* pixel = (SColor*)&image.m_pixels[y * image.m_pitch];
        for (size_t x = 0; x < image.m_width; ++x)
        {
            lambda(*pixel, pixelIndex, x, y);
            ++pixel;
            ++pixelIndex;
        }
    }
}

//======================================================================================
template <typename LAMBDA>
void ImageForEachPixel (const SImageData& image, const LAMBDA& lambda)
{
    size_t pixelIndex = 0;
    for (size_t y = 0; y < image.m_height; ++y)
    {
        SColor* pixel = (SColor*)&image.m_pixels[y * image.m_pitch];
        for (size_t x = 0; x < image.m_width; ++x)
        {
            lambda(*pixel, pixelIndex, x, y);
            ++pixel;
            ++pixelIndex;
        }
    }
}

//======================================================================================
template <typename LAMBDA>
void ImageForBlock (const SImageData& image, size_t blockSize, size_t x, size_t y, const LAMBDA& lambda)
{
    size_t pixelIndex = 0;
    for (size_t iy = 0; iy < blockSize; ++iy)
    {
        SColor* pixel = (SColor*)&image.m_pixels[(iy+iy) * image.m_pitch + x * 3];
        for (size_t ix = 0; ix < blockSize; ++ix)
        {
            lambda(*pixel, pixelIndex);
            ++pixel;
            ++pixelIndex;
        }
    }
}

bool ImageSave (const SImageData &image, const char *fileName);
bool ImageLoad (const char *fileName, SImageData& imageData);
void ImageInit (SImageData& image, size_t width, size_t height);

void ImageConvertToLuma (SImageData& image);
void ImageCombine2 (const SImageData& imageA, const SImageData& imageB, SImageData& result);
void ImageCombine3 (const SImageData& imageA, const SImageData& imageB, const SImageData& imageC, SImageData& result);
void ImageCombine5 (const SImageData& imageA, const SImageData& imageB, const SImageData& imageC, const SImageData& imageD, const SImageData& imageE, SImageData& result);
void ImageDither (const SImageData& sourceImage, const SImageData& noiseImage, SImageData& result);  // assumes sourceImage greyscale (aka is converted to luma u8)

void ImageGaussianBlur (const SImageData& sourceImage, SImageData& destImage, float blurSigma);

void ImageDFT (const SImageData &srcImage, SImageDataComplex &destImage);
void GetMagnitudeData (const SImageDataComplex& srcImage, SImageData& destImage);
