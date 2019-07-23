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

    SDL_LockSurface(image);
    uint32_t *pos = (uint32_t *)pixel_pos(image, x, y);

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
    std::cout << "Load texture: " << name << std::endl;

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
    SDL_FreeSurface(image);
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

    return pixels_[y * width_ + x];
}
