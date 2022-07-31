#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(p0.x-p1.x)<std::abs(p0.y-p1.y)) {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }
    if (p0.x>p1.x) {
        std::swap(p0, p1);
    }

    for (int x=p0.x; x<=p1.x; x++) {
        float t = (x-p0.x)/(float)(p1.x-p0.x);
        int y = p0.y*(1.-t) + p1.y*t;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

Vec3f World2Screen(Vec3f pts){
    return Vec3f(int((pts.x+1.)*width/2+.5), int((pts.y+1.)*height/2+.5), pts.z);
}

bool intriangle(Vec2i *t, Vec2f p){
    Vec2i AB(t[1].x - t[0].x, t[1].y - t[0].y);
    Vec2i BC(t[2].x - t[1].x, t[2].y - t[1].y);
    Vec2i CA(t[0].x - t[2].x, t[0].y - t[2].y);

    Vec2i AP(p.x - t[0].x, p.y - t[0].y);
    Vec2i BP(p.x - t[1].x, p.y - t[1].y);
    Vec2i CP(p.x - t[2].x, p.y - t[2].y);

    Vec3i CrossProduct = Vec3i(AB^AP, BC^BP, CA^CP);;
    return (CrossProduct.x >= 0 && CrossProduct.y >= 0 && CrossProduct.z >= 0);
}



void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
    int min_x, max_x, min_y, max_y;
    min_x = std::min(std::min(t0.x, t1.x), t2.x);
    max_x = std::max(std::max(t0.x, t1.x), t2.x);
    min_y = std::min(std::min(t0.y, t1.y), t2.y);
    max_y = std::max(std::max(t0.y, t1.y), t2.y);
    Vec2i t[3] = {t0, t1, t2};
    for(int i = min_x; i <= max_x; ++i) {
        for(int j = min_y; j <= max_y; ++j) {
            Vec2f p(i, j);
            if(intriangle(t, p)) {
                image.set(i, j, color);
            }
        }
    }
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    TGAImage image(width, height, TGAImage::RGB);
    Vec3f light_dir(0,0,-1);
    for (int i=0; i<model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3];
        Vec3f world_coords[3];

        for (int j=0; j<3; j++) {
            Vec3f v = model->vert(face[j]);
            Vec3f s_v = World2Screen(v);
            screen_coords[j] = Vec2i(s_v.x, s_v.y);
            world_coords[j]  = v;
        }

        Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
        n.normalize();
        float intensity = n*light_dir;
        if (intensity>0) {
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity*255, intensity*255, intensity*255, 255));
        }
    }

    image.flip_vertically(); 
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}
