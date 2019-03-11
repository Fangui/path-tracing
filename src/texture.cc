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

    int height = image->h;
    int width  = image->w;

    this->pixels_.resize(height * width);

    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < width; y++)
        {
            this->set_color(x, y, get_pixel(image, x, y));
        }
    }
}
