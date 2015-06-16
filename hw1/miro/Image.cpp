#include "Miro.h"
#include "Image.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef WIN32
// disable useless warnings
#pragma warning(disable:4996)
#endif

Image * g_image = 0;

Image::Image()
{
    m_pixels = 0;
    m_width = 1;
    m_height = 1;
}

Image::~Image()
{
    if (m_pixels)
        delete [] m_pixels;
}

void Image::resize(int width, int height)
{
    if (m_pixels)
        delete [] m_pixels;
    m_pixels = 0;
    m_pixels = new Pixel[width*height];
    memset(m_pixels, 0, width*height*sizeof(Pixel));
    m_width = width;
    m_height = height;
}

void Image::clear(const Vector3& c)
{
    // should be bg color
    for (int y=0; y<m_height; y++)
        for (int x=0; x<m_width; x++)
            setPixel(x, y, c);
}

// map floating point values to byte values for pixels
unsigned char Map(float r)
{
    float rMap = 255*r;
    unsigned char c = rMap>255?255:(unsigned char)rMap;
    return c;
}

void Image::setPixel(int x, int y, const Vector3& p)
{
    // do some tone mapping
    if (x >= 0 && x < m_width && y < m_height && y >= 0)
    {
        m_pixels[y*m_width+x].r = p.x;
        m_pixels[y*m_width+x].g = p.y;
        m_pixels[y*m_width+x].b = p.z;
    }
}

void Image::setPixel(int x, int y, const Pixel& p)
{
    // do some tone mapping
    if (x >= 0 && x < m_width && y < m_height && y >= 0)
    {
        m_pixels[y*m_width+x]= p;
    }
}

Vector3 Image::getPixel(int x, int y) const
{
    if (x >= 0 && x < m_width && y < m_height && y >= 0)
    {
        return Vector3(m_pixels[y*m_width + x].r, m_pixels[y*m_width + x].g, m_pixels[y*m_width + x].b);
    }
    return Vector3(0.0, 0.0, 0.0);
}

void Image::gammaPixels(float * gammaData)
{
    for (int y = 0; y < m_height; y++)
    {
        for (int x = 0; x < m_width; x++)
        {
            float gammaRed = pow(m_pixels[y*m_width+x].r,1/2.2);
            float gammaGreen = pow(m_pixels[y*m_width+x].g,1/2.2);
            float gammaBlue = pow(m_pixels[y*m_width+x].b,1/2.2);
            gammaData[3 * (m_width*y + x) + 0] = gammaRed;
            gammaData[3 * (m_width*y + x) + 1] = gammaGreen;
            gammaData[3 * (m_width*y + x) + 2] = gammaBlue;
        }
    }
}

void Image::gammaPixels(float * gammaData, int y)
{
    for (int x = 0; x < m_width; x++)
    {
        gammaData[3*(m_width*y+x)+0] = pow(m_pixels[y*m_width+x].r,1/2.2);
        gammaData[3*(m_width*y+x)+1] = pow(m_pixels[y*m_width+x].g,1/2.2);
        gammaData[3*(m_width*y+x)+2] = pow(m_pixels[y*m_width+x].b,1/2.2);
    }
}

void Image::drawScanline(int y)
{
    float * gamma_image = new float[3 * m_width * m_height];
    gammaPixels(gamma_image, y);
    glRasterPos2f(-1, -1 + 2*y / (float)m_height);
    glDrawPixels(m_width, 1, GL_RGB, GL_FLOAT, &gamma_image[3*y*m_width]);
    delete[] gamma_image;
}

void Image::draw()
{
    float * gamma_image = new float[3 * m_width * m_height];
    gammaPixels(gamma_image);
    glDrawPixels(m_width, m_height, GL_RGB, GL_FLOAT, &gamma_image[0]);
    delete[] gamma_image;
}

void Image::writePPM(char* pcFile)
{
    writePPM(pcFile, (float*)m_pixels, m_width, m_height);
}

void Image::writePPM(char *pcFile, float *data, int width, int height)
{
    FILE *fp = fopen(pcFile, "wb");
    if (!fp)
        fprintf(stderr, "Couldn't open PPM file %s for writing\n", pcFile);
    else
    {
        fprintf(fp, "P6\n");
        fprintf(fp, "%d %d\n", width, height );
        fprintf(fp, "255\n" );

        unsigned char * tmp_data = new unsigned char[height * width * 3];
        float * gamma_image = new float[3 * m_width * m_height];
        gammaPixels(gamma_image);
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                tmp_data[3 * (y * width + x) + 0] = Map(gamma_image[3*(y*width + x)+0]);
                tmp_data[3 * (y * width + x) + 1] = Map(gamma_image[3*(y*width + x)+1]);
                tmp_data[3 * (y * width + x) + 2] = Map(gamma_image[3*(y*width + x)+2]);
            }
        }
        delete[] gamma_image;

        // invert image
        int stride = width*3;
        for (int i = height-1; i >= 0; i--)
            fwrite(&tmp_data[stride*i], stride, 1, fp);
        delete[] tmp_data;
        fclose(fp);
    }
}
