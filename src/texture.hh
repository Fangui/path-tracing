#pragma once

#include <vector>
#include <string>

#include "vector.hh"

class Texture
{
    public:
        explicit Texture(const std::string &file);

        int get_width() const
        {
            return this->width_;
        }

        int get_height() const
        {
            return this->height_;
        }

        Vector get_color(int x, int y)
        {
            return this->pixels_[x * this->height_ + y];
        }

        Vector& set_color(int x, int y, Vector c)
        {
            return this->pixels_[x * this->height_ + y] = c;
        }

    private:
        int width_;
        int height_;

        std::vector<Vector> pixels_;
};
