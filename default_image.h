// Simple in-memory image implementation.
#ifndef DEFAULT_IMAGE_H
#define DEFAULT_IMAGE_H

#include "image.h"

namespace sh {

class DefaultImage : public Image {
 public:
  DefaultImage(int width, int height);
  
  int width() const override;
  int height() const override;

  Eigen::Array3f GetPixel(int x, int y) const override;
  void SetPixel(int x, int y, const Eigen::Array3f& v) override;

 private:
  const int width_;
  const int height_;

  std::unique_ptr<Eigen::Array3f[]> pixels_;
};

}  // namespace sh

#endif  // DEFAULT_IMAGE_H
