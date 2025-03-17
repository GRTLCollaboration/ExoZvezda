/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEIGHTFUNCTION_HPP_
#define WEIGHTFUNCTION_HPP_

/*
 * This class defined the weight function needed for unequal-mass fix
 */

class WeightFunction
{
  public:
    WeightFunction() {}

    // Weight function used in https://arxiv.org/abs/2212.08023
    double profile_chi(double coord_x, double coord_y, double coord_z,
                       double radius_width)
    {
        double denom = sqrt(pow(radius_width, 2) + pow(coord_x, 2) +
                            pow(coord_y, 2) + pow(coord_z, 2));
        return 1. / denom;
    }

    // Alternative weight function *UNUSED*
    double compute_weight(double scaledr, int n) const
    {
        double weightfunc;

        if (scaledr <= 1.0)
        {
            if (n == 1)
            {
                weightfunc = 6.0 * ((1.0 / 2.0) * pow((1 - scaledr), 2) -
                                    (1.0 / 3.0) * pow((1 - scaledr), 3));
            }
            if (n == 2)
            {
                weightfunc = 30.0 * ((1.0 / 3.0) * pow((1 - scaledr), 3) -
                                     (1.0 / 2.0) * pow((1 - scaledr), 4) +
                                     (1.0 / 5.0) * pow((1 - scaledr), 5));
            }
            if (n == 3)
            {
                weightfunc = 140.0 * ((1.0 / 4.0) * pow((1 - scaledr), 4) -
                                      (3.0 / 5.0) * pow((1 - scaledr), 5) +
                                      (1.0 / 2.0) * pow((1 - scaledr), 6) -
                                      (1.0) / (7.0) * pow((1 - scaledr), 7));
            }

            if (n == 4)
            {
                weightfunc = 630.0 * ((1.0 / 5.0) * pow((1 - scaledr), 5) -
                                      (2.0 / 3.0) * pow((1 - scaledr), 6) +
                                      (6.0 / 7.0) * pow((1 - scaledr), 7) -
                                      (1.0) / (2.0) * pow((1 - scaledr), 8) +
                                      (1.0) / (9.0) * pow((1 - scaledr), 9));
            }

            if (n == 5)
            {
                weightfunc = 2772.0 * ((1.0 / 6.0) * pow((1 - scaledr), 6) -
                                       (5.0 / 7.0) * pow((1 - scaledr), 7) +
                                       (5.0 / 4.0) * pow((1 - scaledr), 8) -
                                       (10.0) / (9.0) * pow((1 - scaledr), 9) +
                                       (1.0) / (2.0) * pow((1 - scaledr), 10) -
                                       (1.0) / (11.0) * pow((1 - scaledr), 11));
            }
        }

        else
        {
            MayDay::Error(
                "You have requested n parameter larger than implemented!");
            weightfunc = 0.0;
        }

        return weightfunc;
    }
};

#endif /* WEIGHTFUNCTION_HPP_ */
