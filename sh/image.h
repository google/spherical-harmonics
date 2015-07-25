// Copyright 2015 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
