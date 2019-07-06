/*

  == Explanation and rational of testing ==
   Tested each functions to ensure they work for basic cases. Simple tests for
   relatively simple functions. As long as the constructors work as intended,
   most functions dont require much testing. Base and edge cases have mostly
   been skipped so all possibilities have not been covered. Cases where
   exceptions or assertions occur are also not tested.

*/
#include <utility>

#include "assignments/ev/euclidean_vector.h"
#include "catch.h"

SCENARIO("Constructors") {
  GIVEN("an integer x") {
    int x = 2;

    WHEN("a EuclideanVector is constructed using x") {
      EuclideanVector ev{x};

      THEN("the EuclideanVector contains x 0.0's") {
        REQUIRE(std::fabs(ev[0] - 0.0) < 2.22045e-016);
        REQUIRE(std::fabs(ev[1] - 0.0) < 2.22045e-016);
      }
    }
  }

  GIVEN("an integer x and a double y") {
    int x = 2;
    double y = 9.2;

    WHEN("a EuclideanVector is constructed") {
      EuclideanVector ev{x, y};

      THEN("the EuclideanVector contains x y's") {
        REQUIRE(std::fabs(ev[0] - 9.2) < 2.22045e-016);
        REQUIRE(std::fabs(ev[1] - 9.2) < 2.22045e-016);
      }
    }
  }

  GIVEN("a vector") {
    std::vector<double> v{2.3, 5.4};

    WHEN("a EuclideanVector is constructed using the vector") {
      EuclideanVector ev{v.begin(), v.end()};

      THEN("the EuclideanVector contains the same elements as the vector") {
        REQUIRE(std::fabs(ev[0] - 2.3) < 2.22045e-016);
        REQUIRE(std::fabs(ev[1] - 5.4) < 2.22045e-016);
      }
    }
  }

  GIVEN("a EuclideanVector") {
    EuclideanVector ev{6, 1.2};

    WHEN("a EuclideanVector is copy constructed using another EuclideanVector") {
      EuclideanVector c{ev};

      THEN("the new EuclideanVector is a copy of the original EuclideanVector") {
        REQUIRE(c == ev);
      }
    }

    WHEN("a EuclideanVector is move constructed using another EuclideanVector") {
      EuclideanVector m{std::move(ev)};

      THEN("the new EuclideanVector contains the values of the moved EuclideanVector") {
        REQUIRE(m == EuclideanVector{6, 1.2});
      }
    }
  }
}

SCENARIO("Assignments") {
  GIVEN("two EuclideanVectors ev1 & ev2") {
    EuclideanVector ev1{3, 9.53}, ev2{2, 2.22};

    WHEN("ev1 is copy assigned to ev2") {
      ev1 = ev2;

      THEN("ev1 contains the same values as ev2") { REQUIRE(ev1 == ev2); }
    }

    WHEN("ev2 is move assigned to ev1") {
      ev2 = std::move(ev1);

      THEN("ev2 contains ev1's elements") { REQUIRE(ev2 == EuclideanVector{3, 9.53}); }
    }
  }
}

SCENARIO("Type Conversions") {
  GIVEN("a EuclideanVector") {
    std::vector<double> v{4.4, 2.1, 9.5};
    EuclideanVector ev{v.begin(), v.end()};

    WHEN("converting the EuclideanVector to a vector") {
      auto res = std::vector<double>{ev};

      THEN("the vector has the same elements") {
        REQUIRE(static_cast<int>(res.size()) == ev.GetNumDimensions());
        REQUIRE(res == v);
      }
    }

    WHEN("converting the EuclideanVector to a list") {
      auto res = std::list<double>{ev};

      THEN("the vector has the same elements") {
        REQUIRE(static_cast<int>(res.size()) == ev.GetNumDimensions());
        REQUIRE(res == std::list<double>{v.begin(), v.end()});
      }
    }
  }
}

SCENARIO("Using Subscript") {
  GIVEN("a EuclideanVector and an integer") {
    std::vector<double> v{4.4, 2.1, 9.5};
    EuclideanVector ev{v.begin(), v.end()};
    int i = 2;

    WHEN("getting the value at the integer using subscript") {
      auto res = ev[i];

      THEN("the value at that index is returned") { REQUIRE(std::fabs(res - 9.5) < 2.22045e-016); }
    }

    WHEN("setting the value at the integer using subscript") {
      ev[i] = 9.99;

      THEN("the value at that index is changed") { REQUIRE(ev[i] == 9.99); }
    }
  }
}

SCENARIO("Using Operators using two EuclideanVectors") {
  GIVEN("two different EuclideanVectors of the same dimension") {
    std::vector<double> v1{3.5, 1.2, 1.8}, v2{2.5, 1.8, 3.2};
    std::vector<double> res_add{3.5 + 2.5, 1.2 + 1.8, 1.8 + 3.2},
        res_sub{3.5 - 2.5, 1.2 - 1.8, 1.8 - 3.2};
    EuclideanVector ev1{v1.begin(), v1.end()}, ev2{v2.begin(), v2.end()};
    EuclideanVector res_ev_add{res_add.begin(), res_add.end()},
        res_ev_sub{res_sub.begin(), res_sub.end()};

    WHEN("checking equality between 2 different vectors") {
      bool b1 = (v1 == v2), b2 = (v1 != v2);

      THEN("the result is false") {
        REQUIRE(!b1);
        REQUIRE(b2);
      }
    }

    WHEN("checking equality between 2 same vectors") {
      bool b1 = (v1 == v1), b2 = (v1 != v1);

      THEN("the result is true") {
        REQUIRE(b1);
        REQUIRE(!b2);
      }
    }

    WHEN("adding the two") {
      auto res = ev1 + ev2;
      ev1 += ev2;

      THEN("the resulting vector is the point-wise sum") {
        REQUIRE(res == res_ev_add);
        REQUIRE(ev1 == res_ev_add);
      }
    }

    WHEN("subtracting the two") {
      auto res = ev1 - ev2;
      ev1 -= ev2;

      THEN("the resulting vector is the point-wise sum") {
        REQUIRE(res == res_ev_sub);
        REQUIRE(ev1 == res_ev_sub);
      }
    }

    WHEN("finding the dot-product of the two") {
      auto res = ev1 * ev2;

      THEN("the resulting vector is the point-wise sum") {
        REQUIRE(std::fabs(3.5 * 2.5 + 1.2 * 1.8 + 1.8 * 3.2 - res) < 2.22045e-016);
      }
    }
  }
}

SCENARIO("Using Operators that use a EuclideanVector and a Scalar") {
  GIVEN("a EuclideanVector and a scalar") {
    int scalar = 10;
    std::vector<double> v{0.5, 6.8, 3.2};
    std::vector<double> res_mul{0.5 * scalar, 6.8 * scalar, 3.2 * scalar},
        res_div{0.5 / scalar, 6.8 / scalar, 3.2 / scalar};
    EuclideanVector ev{v.begin(), v.end()};
    EuclideanVector res_ev_mul{res_mul.begin(), res_mul.end()},
        res_ev_div{res_div.begin(), res_div.end()};

    WHEN("multiplying the two") {
      auto res = ev * scalar;
      ev *= scalar;

      THEN("each element in the resulting vector is scaled by scalar") {
        REQUIRE(res == res_ev_mul);
        REQUIRE(ev == res_ev_mul);
      }
    }

    WHEN("dividing the two") {
      auto res = ev / scalar;
      ev /= scalar;

      THEN("each element in the resulting vector is scaled by 1/scalar") {
        REQUIRE(res == res_ev_div);
        REQUIRE(ev == res_ev_div);
      }
    }
  }
}

SCENARIO("EuclideanVector Methods") {
  GIVEN("a EuclideanVector with different values") {
    std::vector<double> v{2.0, 6.0, 3.0};
    EuclideanVector ev{v.begin(), v.end()};

    WHEN("getting its dimension") {
      auto dim = ev.GetNumDimensions();

      THEN("the dimension is returned") { REQUIRE(dim == 3); }
    }

    GIVEN("one of its index") {
      int index = 1;

      WHEN("getting the value at that index") {
        auto value = ev.at(index);

        THEN("the value at the index is returned") {
          REQUIRE(std::fabs(value - 6.0) < 2.22045e-016);
        }
      }

      WHEN("setting the value at that index") {
        ev.at(index) = 2.0;

        THEN("the value at the index is returned") {
          REQUIRE(std::fabs(ev.at(index) - 2.0) < 2.22045e-016);
        }
      }
    }

    WHEN("getting the Euclidean Norm") {
      auto norm = ev.GetEuclideanNorm();

      THEN("the euclidean norm is returned") { REQUIRE(std::fabs(norm - 7.0) < 2.22045e-016); }
    }

    WHEN("getting its unit vector") {
      auto uv = ev.CreateUnitVector();
      std::vector<double> uv_v{2.0 / 7.0, 6.0 / 7.0, 3.0 / 7.0};
      EuclideanVector correct_uv{uv_v.begin(), uv_v.end()};

      THEN("the unit vector is returned") { REQUIRE(uv == correct_uv); }
    }
  }
}
