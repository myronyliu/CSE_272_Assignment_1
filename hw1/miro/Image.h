#ifndef IMAGE_H_INCLUDED
#define IMAGE_H_INCLUDED

#include "Vector3.h"

class Image
{
public:
    struct Pixel
    {
        float r, g, b;
        Pixel(float ir, float ig, float ib) {set(ir, ig, ib);}
        Pixel() : r(0), g(0), b(0) {}
        void set(float ir, float ig, float ib) {r = ir; g = ig; b = ib;}
    };

    Image();
    ~Image();

    void resize(int width, int height);
    void setPixel(int x, int y, const Vector3& p);
    void setPixel(int x, int y, const Pixel& p);

    Vector3 getPixel(int x, int y) const;


    void Image::gammaPixels(float * gammaData);
    void Image::gammaPixels(float * gammaData, int y);
    void draw();
    void drawScanline(int y);
    void clear(const Vector3& c);
    void writePPM(char* pcFile); // write data to a ppm image file
    void writePPM(char *pcName, float *data, int width, int height);

    float* getCharPixels()  {return (float*)m_pixels;}
    int width() const               {return m_width;}
    int height() const              {return m_height;}

private:
    Pixel* m_pixels;
    int m_width;
    int m_height;
};

extern Image * g_image;

#endif // IMAGE_H_INCLUDED
