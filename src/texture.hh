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

        unsigned get_size() const
        {
            return pixels_.size();
        }

        Vector get_color(int x, int y) const
        {
            return this->pixels_[x * this->height_ + y];
        }

        Vector get_color(double u, double v) const;


        void set_color(int x, int y, Vector c)
        {
            this->pixels_[x * this->height_ + y] = c;
        }

        int width_;
        int height_;

        std::vector<Vector> pixels_;
};
