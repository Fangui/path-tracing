#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "vector.hh"

std::vector<Vector> load_vector(std::string &path)
{
    std::vector <Vector> vect;
    std::ifstream in(path);
    if (!in.is_open())
    {
        std::cerr << "Path not found: " << path << '\n';
        exit(2);
    }

    std::string line;
    std::getline(in, line); // remve hearder
    std::getline(in, line); // width height
    std::getline(in, line); // max

    size_t height = 512;
    double a;
    double b;
    double c;

    while (std::getline(in, line))
    {
        std::istringstream strin(line);
        for (size_t i = 0; i < height; ++i)
        {

            strin >> a;
            strin >> b;
            strin >> c;

            a /= 255;
            b /= 255;
            c /= 255;
            vect.push_back(Vector(a, b, c));
        }
    }

    std::cout << vect.size() << std::endl;
    return vect;
}

int write_ppm(const std::string &out_path, const std::vector<Vector> &vect,
        int width, int height)
{
    std::ofstream out (out_path);
    unsigned index = 0;
    if (out.is_open())
    {
        out << "P3\n";
        out << width << " " << height << '\n';
        out << 255 << '\n';

        for (int i = 0; i < width; ++i)
        {
            for (int j = 0; j < height; ++j)
            {
                int r = vect[index][0] * 255.0;
                int g = vect[index][1] * 255.0;
                int b = vect[index++][2] * 255.0;

                r = std::min(r, 255);
                g = std::min(g, 255);
                b = std::min(b, 255);
                out << r << " " << g << " " << b << "  ";
            }
            out << '\n';
        }
        std::cout << "Create " + out_path + " file\n";
    }
    else
    {
        std::cerr << "Error while write in " << out_path << '\n';
        return 1;
    }
    return 0;
}

std::vector<Vector> merge(std::vector<Vector> vect1, 
                          std::vector<Vector> vect2, double coef)
{
    for (size_t i = 0; i < vect1.size(); ++i)
        vect1[i] = vect1[i] * coef + vect2[i] * (1 - coef);

    return vect1;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " path1 path2 coef\n";
        return 1;
    }

    std::string path1 = argv[1];
    std::string path2 = argv[2];
    double coef = 0.5;
    if (argc > 3)
        coef = atof(argv[3]);

    if (coef > 1)
    {
        std::cerr << "Error, coef must be between 0 and 1\n";
        return 2;
    }

    std::cout <<  "coef: " << coef << std::endl;

    auto vect1 = load_vector(path1);
    auto vect2 = load_vector(path2);

    if (vect1.size() != vect2.size())
    {
        std::cerr << "Error vector must has same size\n";
        return 2;
    }

    vect1 = merge(vect1, vect2, coef);
    return write_ppm("sum.ppm", vect1, 512, 512);
}
