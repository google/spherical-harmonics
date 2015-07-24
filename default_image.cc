#include "default_image.h"

namespace sh {

DefaultImage::DefaultImage(int width, int height) : width_(width), height_(height) {
  pixels_.reset(new Eigen::Array3f[width * height]);
}

int DefaultImage::width() const { return width_; }

int DefaultImage::height() const { return height_; }

Eigen::Array3f DefaultImage::GetPixel(int x, int y) const {
  int index = x + y * width_;
  return pixels_[index];
}

void DefaultImage::SetPixel(int x, int y, const Eigen::Array3f& v) {
  int index = x + y * width_;
  pixels_[index] = v;
}

}  // namespace sh
