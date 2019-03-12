#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#include "texture.hh"

void *pixel_pos(SDL_Surface *image, int x, int y)
{
    uint8_t *p = (uint8_t *) image->pixels
                        + y * image->pitch
                        + x * image->format->BytesPerPixel;

    return (void *) p;
}

Vector get_pixel(SDL_Surface *image, int x, int y)
{
    uint8_t r = 0.;
    uint8_t g = 0.;
    uint8_t b = 0.;

    uint32_t *pos = (uint32_t *) pixel_pos(image, x, y);

    SDL_GetRGB(*pos, image->format, &r, &g, &b);
    return Vector{r / 255.0f, g / 255.0f, b / 255.0f};
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
#include <iostream>
int main(int argc, char *argv[])
{
    if (argc < 2)
        return 1;

    Texture t(argv[1]);
    std::cout << t.get_height() <<std::endl;
    std::cout << t.get_width() <<std::endl;
    for (int i = 0; i < t.get_width(); ++i)
    {
        for (int j = 0; j < t.get_height(); ++j)
        {
            std::cout << t.get_color(i, j) << std::endl;
        }
    }
    return 0;
}*/
