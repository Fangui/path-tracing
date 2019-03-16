#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#include "texture.hh"

/*
Uint32 pixel_pos(SDL_Surface *surface, int x, int y)
{
    int bpp = surface->format->BytesPerPixel;
    Uint8 *p = (Uint8 *)surface->pixels + y * surface->pitch + x * bpp;

    switch(bpp) {
        case 1:
            return *p;

        case 2:
            return *(Uint16 *)p;

        case 3:
            if(SDL_BYTEORDER == SDL_BIG_ENDIAN)
                return p[0] << 16 | p[1] << 8 | p[2];
            else
                return p[0] | p[1] << 8 | p[2] << 16;

        case 4:
            return *(Uint32 *)p;

        default:
}*/

void *pixel_pos(SDL_Surface *image, int x, int y)
{
    uint8_t *p = (uint8_t *) image->pixels
                        + (x * image->h
                                                                                        + y) * image->format->BytesPerPixel;

    return (void *) p;
}

Vector get_pixel(SDL_Surface *image, int x, int y)
{
    uint8_t r = 0.;
    uint8_t g = 0.;
    uint8_t b = 0.;

    SDL_LockSurface(image);
    uint32_t *pos = (uint32_t*)pixel_pos(image, x, y);

    SDL_UnlockSurface(image); 
    SDL_GetRGB(*pos, image->format, &r, &g, &b);
    return Vector((double)r / 255., (double)g / 255., (double)b / 255.);
}

Texture::Texture(const std::string &name)
{
    SDL_Surface *image = IMG_Load(name.c_str());

    if (!image)
    {
        throw std::runtime_error("cannot open image " + name);
    }

    this->height_ = image->h;
    this->width_  = image->w;

    this->pixels_.resize(height_ * width_);

    for (int x = 0; x < width_; ++x)
    {
        for (int y = 0; y < height_; ++y)
        {
            this->set_color(x, y, get_pixel(image, x, y));
        }
    }
}

Vector Texture::get_color(double u, double v) const
{
    u = fabs(u);
    v = fabs(v);

    while (u > 1)
        u -= 1;
    while (v > 1)
        v -= 1;

    int x = u * (width_ - 1);
    int y = (1 - v) * (height_ - 1);
    
    return pixels_[x * height_ + y];
/*

    int m_Width = this->width_;
    int m_Height = this->height_;

    double fu = (a_U + 1000.5f) * m_Width;
    double fv = (a_V + 1000.0f) * m_Width;
    int u1 = ((int)fu) % m_Width;
    int v1 = ((int)fv) % m_Height;
    int u2 = (u1 + 1) % m_Width;
    int v2 = (v1 + 1) % m_Height;
    double fracu = fu - floorf( fu );
    double fracv = fv - floorf( fv );
    // calculate weight factors
    double w1 = (1 - fracu) * (1 - fracv);
    double w2 = fracu * (1 - fracv);
    double w3 = (1 - fracu) * fracv;
    double w4 = fracu *  fracv;
    // fetch four texels
    Vector c1 = pixels_[u1 + v1 * m_Width];
    Vector c2 = pixels_[u2 + v1 * m_Width];
    Vector c3 = pixels_[u1 + v2 * m_Width];
    Vector c4 = pixels_[u2 + v2 * m_Width];
    // scale and sum the four colors
    return c1 * w1 + c2 * w2 + c3 * w3 + c4 * w4;
*/

}



/*
#include <fstream>
#include <iostream>
int write_ppm(const std::string &out_path, const Texture &t,
                int width, int height)
{
    std::ofstream out (out_path);
    unsigned index = 0;
    if (out.is_open())
    {
            out << "P3\n";
            out << width << " " << height << '\n';
            out << 255 << '\n';
    
            for (double i = 0; i < 1; i += 1 / (double)width)
            {
                        for (double j = 0; j < 1; j += 1 / (double)height)
                        {
                            Vector color = t.get_color(i, j);
                            int r = color[0] * 255.0;
                            int g = color[1] * 255.0;
                            int b = color[2] * 255.0;
                            out << r << " " << g << " " << b << "  ";
                                    }
                        out << '\n';
                    }
        }
    else
    {
            std::cerr << "Error while write in " << out_path << '\n';
            return 1;
        }
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
        return 1;

    Texture t(argv[1]);
    std::cout << t.get_height() <<std::endl;
    std::cout << t.get_width() <<std::endl;
    
    std::cout << t.get_size() << std::endl;
    write_ppm("test", t, t.get_width(), t.get_height());
    return 0;
}*/
