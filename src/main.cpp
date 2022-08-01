#include <vector>
#include <cmath>
#include <limits>
#include <tuple>
#include <memory>

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"


const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
std::unique_ptr<Model> model;
const int width = 800;
const int height = 800;

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y)) {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }
    if (p0.x > p1.x) {
        std::swap(p0, p1);
    }

    for (int x = p0.x; x <= p1.x; x++) {
        float t = (x - p0.x) / (float)(p1.x - p0.x);
        int y = p0.y * (1. - t) + p1.y * t;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

Vec3f World2Screen(Vec3f pts) {
    return Vec3f(int((pts.x + 1.) * width / 2 + .5), int((pts.y + 1.) * height / 2 + .5), pts.z);
}

bool intriangle(Vec2i *t, Vec2f p) {
    Vec2i AB(t[1].x - t[0].x, t[1].y - t[0].y);
    Vec2i BC(t[2].x - t[1].x, t[2].y - t[1].y);
    Vec2i CA(t[0].x - t[2].x, t[0].y - t[2].y);

    Vec2i AP(p.x - t[0].x, p.y - t[0].y);
    Vec2i BP(p.x - t[1].x, p.y - t[1].y);
    Vec2i CP(p.x - t[2].x, p.y - t[2].y);

    Vec3i CrossProduct = Vec3i(AB ^ AP, BC ^ BP, CA ^ CP);;
    return (CrossProduct.x >= 0 && CrossProduct.y >= 0 && CrossProduct.z >= 0);
}

std::tuple <float, float, float> computeBarycentric(Vec3f *triangle, int x, int y) {
    float alpha = ((x - triangle[2].x) * (triangle[0].y - triangle[1].y) - (triangle[0].x - triangle[1].x) * (y - triangle[2].y)) / ((triangle[0].x - triangle[2].x) * (triangle[0].y - triangle[1].y) - (triangle[0].x - triangle[1].x) * (triangle[0].y - triangle[2].y));
    float beta = ((x - triangle[2].x) - alpha * (triangle[0].x - triangle[2].x)) / (triangle[0].x - triangle[1].x);
    float gamma = 1 - alpha - beta;
    return {alpha, beta, gamma};
}

void triangle(Vec2i *s_coords, Vec3f *w_coords, TGAImage &image, TGAColor color, float *zbuffer) {
    int min_x, max_x, min_y, max_y;

    // min_x = std::min({s_coords[0].x, s_coords[1].x, s_coords[2].x});
    min_x = std::min(std::min(s_coords[0].x, s_coords[1].x), s_coords[2].x);
    max_x = std::max(std::max(s_coords[0].x, s_coords[1].x), s_coords[2].x);
    min_y = std::min(std::min(s_coords[0].y, s_coords[1].y), s_coords[2].y);
    max_y = std::max(std::max(s_coords[0].y, s_coords[1].y), s_coords[2].y);
    for (int i = min_x; i <= max_x; ++i) {
        for (int j = min_y; j <= max_y; ++j) {
            Vec2f p(i, j);
            if (intriangle(s_coords, p)) {
                float w_i = (i - .5) * 2 / width - 1., w_j = (j - .5) * 2 / width - 1.;
                auto [alpha, beta, gamma] = computeBarycentric(w_coords, w_i + .5, w_j + .5);
                float z = alpha * w_coords[0].z + beta* w_coords[1].z + gamma * w_coords[2].z;
                if (z > zbuffer[i + j * width]) {
                    zbuffer[i + j * width] = z;
                    image.set(i, j, color);
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    if (2 == argc) {
        model = std::make_unique<Model>(argv[1]);
    } else {
        model = std::make_unique<Model>("obj/african_head.obj");
    }

    TGAImage image(width, height, TGAImage::RGB);
    Vec3f light_dir(0, 0, -1);

    //zbuffer init
    auto zbuffer = std::make_unique<float[]>(width * height);

    for (int i = 0; i < width * height; ++i) {
        zbuffer[i] = -std::numeric_limits<float>::max();
    }

    //triangle
    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3];
        Vec3f world_coords[3];

        for (int j = 0; j < 3; j++) {
            Vec3f v = model->vert(face[j]);
            Vec3f s_v = World2Screen(v);
            screen_coords[j] = Vec2i(s_v.x, s_v.y);
            world_coords[j] = v;
        }

        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        n.normalize();
        float intensity = n * light_dir;
        if (intensity > 0) {
            triangle(screen_coords, world_coords, image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255), zbuffer.get());
        }
    }

    image.flip_vertically();
    image.write_tga_file("output.tga");
    return 0;
}
