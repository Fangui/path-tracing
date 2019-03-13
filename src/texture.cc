#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#include "texture.hh"

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

    uint32_t *pos = (uint32_t *) pixel_pos(image, x, y);

    SDL_GetRGB(*pos, image->format, &r, &g, &b);
    return Vector(r / 255., g / 255., b / 255.);
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
