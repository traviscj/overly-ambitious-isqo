// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)
// test.cpp
#include <UnitTest++.h>

TEST(FailSpectacularly)
{
  CHECK(false);
}

int main()
{
  return UnitTest::RunAllTests();
}
