//
// Created by LEI XU on 5/12/19.
//

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>

#include "global.hpp"
#include "Object.hpp"
#include "Vector.hpp"
#include "Sphere.hpp"
#include "Triangle.hpp"
#include "Scene.hpp"
#include "Light.hpp"
#include "AreaLight.hpp"

inline float deg2rad(const float &deg)
{ return deg * M_PI/180.0; }

const float  EPSILON = 0.00001;


// The main render function. This where we iterate over all pixels in the image, generate
// primary rays and cast these rays into the scene. The content of the framebuffer is
// saved to a file.
void render(Scene& scene)
{
    Vector3f *framebuffer = new Vector3f[scene.width * scene.height];
    Vector3f *pix = framebuffer;
    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f orig(-1, 5, 10);
    for (uint32_t j = 0; j < scene.height; ++j) {
        for (uint32_t i = 0; i < scene.width; ++i) {
            float x = (2 * ((i + 0.5) / scene.width) - 1) * scale * imageAspectRatio;
            float y = (1 - 2 * ((j + 0.5) / scene.height)) * scale;

            Vector3f dir = Vector3f(x, y, -1);
            normalize(dir);
            *(pix++) = scene.castRay(Ray(orig, dir), 0);
        }
        UpdateProgress(j / (float)scene.height);
    }

    // save framebuffer to file
    std::ofstream ofs;
    ofs.open("./out.ppm");
    ofs << "P6\n" << scene.width << " " << scene.height << "\n255\n";
    for (uint32_t i = 0; i < scene.height * scene.width; ++i) {
        char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
        char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
        char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
        ofs << r << g << b;
    }

    ofs.close();

    delete [] framebuffer;
}

// In the main function of the program, we create the scene (create objects and lights)
// as well as set the options for the render (image width and height, maximum recursion
// depth, field-of-view, etc.). We then call the render function().
int main(int argc, char **argv)
{
    time_t start, stop;

    Scene scene;
    scene.width = 1280;
    scene.height = 960;
    scene.fov = 90;
    scene.backgroundColor = Vector3f(0.235294, 0.67451, 0.843137);
    scene.maxDepth = 5;

    MeshTriangle bunny("../models/bunny/bunny.obj");

    scene.Add(&bunny);
    scene.Add(std::make_unique<Light>(Vector3f(-20, 70, 20), 1));
    scene.Add(std::make_unique<Light>(Vector3f(20, 70, 20), 1));
    scene.buildBVH();

    //render
    time(&start);
    render(scene);
    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff/3600;
    int mins = ((int)diff/60)-(hrs*60);
    int secs = (int)diff-(hrs*3600)-(mins*60);

    printf("\rRender complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n", hrs, mins, secs);

    return 0;
}
