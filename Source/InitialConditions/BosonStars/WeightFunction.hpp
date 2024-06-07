#ifndef WEIGHTFUNCTION_HPP_
#define WEIGHTFUNCTION_HPP_

class WeightFunction
{
    public:

    WeightFunction(){}

    double compute_weight(double scaledr, int n) const
    {   
        double weightfunc;
      
        if (scaledr <= 1.0)
        {   
            if (n==1)
            {
                weightfunc = 6.0 * ((1.0/2.0) * pow((1-scaledr), 2) - (1.0/3.0) * pow((1-scaledr), 3));
            }
            if (n==2)
            {
                weightfunc = 30.0 * ((1.0/3.0) * pow((1-scaledr), 3) - (1.0/2.0) * pow((1-scaledr), 4) + (1.0/5.0) * pow((1-scaledr), 5));
            }
            if (n==3)
            {
                weightfunc = 140.0 * ( (1.0/4.0)*pow((1-scaledr),4) - (3.0/5.0)*pow((1-scaledr),5) + (1.0/2.0)* pow((1-scaledr),6) - (1.0)/(7.0)*pow((1-scaledr),7));
            }

            if (n==4)
            {
                weightfunc = 630.0 * ( (1.0/5.0)*pow((1-scaledr),5) - (2.0/3.0)*pow((1-scaledr),6) + (6.0/7.0)* pow((1-scaledr),7) - (1.0)/(2.0)*pow((1-scaledr),8) + (1.0)/(9.0)*pow((1-scaledr),9));
            }

             if (n==5)
            {
                weightfunc = 2772.0 * ( (1.0/6.0)*pow((1-scaledr),6) - (5.0/7.0)*pow((1-scaledr),7) + (5.0/4.0)* pow((1-scaledr),8) - (10.0)/(9.0)*pow((1-scaledr),9) + (1.0)/(2.0)*pow((1-scaledr),10) - (1.0)/(11.0)*pow((1-scaledr),11));
            }
            
	        // weightfunc = 4808643120 * ((1.0/16.0) * pow(1-scaledr, 16) - (15.0/17.0) * pow(1-scaledr, 17) + (35.0/6.0) * pow(1-scaledr, 18) - (455.0/19.0) * pow(1-scaledr, 19) + (273.0/4.0) * pow(1-scaledr, 20) - 143.0 * pow(1-scaledr, 21) + (455.0/2.0) * pow(1-scaledr, 22) - (6435.0/23.0) * pow(1-scaledr, 23) + (2145.0/8.0) * pow(1-scaledr, 24) - (1001.0/5.0) * pow(1-scaledr, 25) + (231.0/2.0) * pow(1-scaledr, 26) - (455.0/9.0) * pow(1-scaledr, 27) + (65.0/4.0) * pow(1-scaledr, 28) - (105.0/29.0) * pow(1-scaledr, 29) + (1.0/2.0) * pow(1-scaledr, 30) - (1.0/31.0) * pow(1-scaledr, 31));  
	}

        else
        {
            weightfunc = 0.0;
        }

    return weightfunc;
    }

    double stretching_factor(double coord_x, double coord_y, double alpha_zero)
    {
        double sin_half_angle_squared = (1.0 / 2.0) * (1.0 - (coord_x) / (sqrt(pow(coord_x, 2) + pow(coord_y, 2))));

        double stretch = alpha_zero + (1.0 - alpha_zero) * sin_half_angle_squared;

        return stretch;
    }

    double stretching_factor2(double coord_x, double coord_y, double alpha_zero)
    {
        double cos_half_angle_squared = (1.0 / 2.0) * (1.0 + (coord_x) / (sqrt(pow(coord_x, 2) + pow(coord_y, 2))));

        double stretch2 = alpha_zero + (1.0 - alpha_zero) * cos_half_angle_squared;

        return stretch2;
    }

    double profile_chi(double coord_x, double coord_y, double coord_z, double radius_width)
    {
        double denom = sqrt(pow(radius_width, 2) + pow(coord_x, 2) + pow(coord_y, 2)+ pow(coord_z, 2));
        return 1. / denom;
    }
};

#endif /* WEIGHTFUNCTION_HPP_ */
