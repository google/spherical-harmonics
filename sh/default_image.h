// Simple in-memory image implementation.
#ifndef SH_DEFAULT_IMAGE_H
#define SH_DEFAULT_IMAGE_H

#include "sh/image.h"

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

#endif  // SH_DEFAULT_IMAGE_H
