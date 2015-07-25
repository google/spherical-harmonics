// A very simple image interface that reports image dimensions and allows
// for setting and getting individual pixels. Implement the interface by
// wrapping the image library of your choice.

#ifndef SH_IMAGE_H
#define SH_IMAGE_H

#include <memory>

#include "Eigen/Dense"

namespace sh {

class Image {
 public:
  Image() {}
  virtual ~Image() {}

  virtual int width() const = 0;
  virtual int height() const = 0;

  virtual Eigen::Array3f GetPixel(int x, int y) const = 0;
  virtual void SetPixel(int x, int y, const Eigen::Array3f& v) = 0;
};

}  // namespace sh

#endif  // IMAGE_H
