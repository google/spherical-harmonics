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
