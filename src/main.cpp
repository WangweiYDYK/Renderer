#include <vector>
#include <cmath>
#include <limits>
#include <tuple>
#include <memory>

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "gl.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);

std::unique_ptr<Model> model;

Vec3f light_dir(1, 1, 1);
Vec3f camera(0, -1, 3);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);

const int width = 800;
const int height = 800;

struct GouraudShader : public IShader{
    Vec3f varying_intensity;

    Vec4f vertex(int iface, int nthvert) override {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;     // transform it to screen coordinates
        varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity
        return gl_Vertex;
    }

    bool fragment(Vec3f bar, TGAColor &color) override {
        float intensity = varying_intensity*bar;
        color = TGAColor(intensity*255, intensity*255, intensity*255, 255);
        return false;
    }
};


int main(int argc, char **argv) {
    if (2 == argc) {
        model = std::make_unique<Model>(argv[1]);
    } else {
        model = std::make_unique<Model>("obj/african_head.obj");
    }

    //matrix init
    light_dir.normalize();
    lookat(camera, center, up);
    viewport(width/8, height/8, width*3/4, height*3/4);
    projection(-1.f /(camera - center).norm());

    TGAImage image(width, height, TGAImage::RGB);
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

    GouraudShader shader;
    for (int i = 0; i < model->nfaces(); i++) {
        Vec4f screen_coords[3];
        for(int j = 0; j < 3; ++j){
            screen_coords[j] = shader.vertex(i, j);
        }
        triangle(screen_coords, shader, image, zbuffer);
    }
    image.flip_vertically();
    zbuffer.flip_vertically();
    zbuffer.write_tga_file("zbuffer.tga");
    image.write_tga_file("output.tga");
    return 0;
}
