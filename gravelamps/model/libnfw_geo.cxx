// Function definitions for calculating the amplification factor for the
// Navarro, Frenk, White (NFW) lens model in the geometric optics apprxoimation
// regime.
//
// Mick Wright 2021

#include "src/nfw.h"

// Function computes the value of the lens equation for a given set of image
// position, source position, and scaling constant of the profile.
//
// Input:
//      double image_position : dimensionless position of the image considered
//      double source_position : dimensionless displacemnt from optical axis
//      double scaling_constant : characteristic scale of the profile
//
// Output:
//      double lens_eq : value of the lens equation with the given parameters
double LensEquation(double image_position,
                    double source_position,
                    double scaling_constant) {
    // Calculate the differential of the lensing potential
    // Setting the stepsize
    const int radix = std::numeric_limits<double>::radix;
    const int mantissa = std::numeric_limits<double>::digits;
    double stepsize = std::pow(1.0/radix, mantissa/2);

    double image_plus = image_position + stepsize;
    double image_minus = image_position - stepsize;

    double potential_plus =\
        std::real(LensingPotential(image_plus, scaling_constant));
    double potential_minus =\
        std::real(LensingPotential(image_minus, scaling_constant));

    double differential = (potential_plus - potential_minus)/(2 * stepsize);

    // Compute the lensing equation
    double lens_eq = image_position - differential - source_position;

    return lens_eq;
}

// Function computes the image positions for a given source position by means
// of root finding the lens equation in the negative, zero, and positive cases.
//
// Input:
//      double source_position : dimensionless displacement from optical axis
//      double scaling_constant : characteristic scale of the profile
//
// Output:
//      std::vector<double> image_positions : vector containing the image
//                                            positions - will either have
//                                            length one or three
std::vector<double> ImagePositions(double source_position,
                                   double scaling_constant) {
    using boost::math::tools::bisect;
    using boost::math::tools::eps_tolerance;

    // Setting desired accuracy to a little less than maximal
    int digits = std::numeric_limits<double>::digits;
    int tol = digits - 3;
    eps_tolerance<double> tolerance(tol);

    // Lambda variant of the Lens Equation to make it single parameter
    auto LensEq =
        [source_position, scaling_constant](double image_position) {
            return LensEquation(image_position,
                                source_position,
                                scaling_constant);
        };

    // Find the three possible roots if they exist. These are a negative root,
    // a near zero root, and a positive root. In the case of the source
    // position being above the critical value, there will be only one
    // root that being the positive one.
    std::vector<double> image_positions;

    std::pair<double, double> root_pos_pair =
        bisect(LensEq, 0.1, 12., tolerance);

    double root_pos = root_pos_pair.first\
                      + (root_pos_pair.second - root_pos_pair.first)/2.;
    image_positions.push_back(root_pos);

    double root_neg = root_pos;
    try {
        std::pair<double, double> root_neg_pair =
            bisect(LensEq, -12., -0.1, tolerance);
        root_neg = root_neg_pair.first
                   + (root_neg_pair.second - root_neg_pair.first)/2.;
    } catch(...) {
    }

    double root_zero = root_pos;
    try {
        std::pair<double, double> root_zero_pair =
            bisect(LensEq, -0.1, 0.1, tolerance);
        root_zero = root_zero_pair.first
                    + (root_zero_pair.second - root_zero_pair.first)/2.;
    } catch(...) {
    }

    if (root_neg != root_pos) {
        image_positions.push_back(root_neg);
    }
    if (root_zero != root_pos) {
        image_positions.push_back(root_zero);
    }

    return image_positions;
}

// Function computes the NFW mass associated with a given image position and
// scaling constant.
//
// Input:
//      double image_position : dimensionless position of the image considered
//      double scaling_constant : characteristic scale of the profile
//
// Output:
//      std::complex<double> image_mass : complex value of the mass associated
//                                        with the image
std::complex<double> ImageMass(double image_position,
                               double scaling_constant) {
    std::complex<double> image_complex = image_position;
    std::complex<double> image_complex_sq = image_position * image_position;

    // Compute the trigonometric term. This has three values depending on the
    // relation of the image position to 1
    std::complex<double> trig_term;

    if (image_position == 1) {
        trig_term = 1.;
    } else if (image_position < 1) {
        std::complex<double> image_position_term = sqrt(1. - image_complex_sq);
        trig_term = atanh(image_position_term)/image_position_term;
    } else {
        std::complex<double> image_position_term = sqrt(image_complex_sq - 1.);
        trig_term = atan(image_position_term)/image_position_term;
    }

    // Compute the logarithmic term
    std::complex<double> log_term = log(image_complex/2.);

    // Compute final result
    std::complex<double> image_mass = scaling_constant * (log_term + trig_term);

    return image_mass;
}

// Function computes the surface density associated with an image at the
// specified position for the given scaling
//
// Input:
//      double image_position : dimensionless position of the image considered
//      double scaling_constant : characteristic scale of the profile
//
// Output:
//      std::complex<double> surface_density : complex value of the surface
//                                             density associated with the
//                                             image position
std::complex<double> SurfaceDensity(double image_position,
                                    double scaling_constant) {
    // Compute the scaling constant term
    double scaling_constant_term = scaling_constant/2.;

    // Compute the position dependent term based on the relation of the image
    // position to one
    std::complex<double> image_complex = image_position;
    std::complex<double> image_complex_sq = image_complex * image_complex;

    std::complex<double> position_dependent_term;
    if ( image_position == 1 ) {
        position_dependent_term = 1./3.;
    } else if ( image_position > 1 ) {
        std::complex<double> position_power_term =\
             -2./pow((image_complex_sq - 1.), 1.5);

        std::complex<double> position_trig_term =\
            atan(sqrt((image_complex - 1.)/(1. + image_complex)));

        std::complex<double> position_frac_term = 1./(image_complex_sq - 1.);

        position_dependent_term = position_power_term * position_trig_term\
                                  + position_frac_term;
    } else {
        std::complex<double> position_power_term =\
            2./pow((1. - image_complex_sq), 1.5);

        std::complex<double> position_trig_term =\
            atanh(sqrt((1. - image_complex)/(image_complex + 1.)));

        std::complex<double> position_frac_term = 1./(1. - image_complex_sq);

        position_dependent_term = position_power_term * position_trig_term\
                                  - position_frac_term;
    }

    // Compute the final surface density value
    std::complex<double> surface_density =\
        scaling_constant_term * position_dependent_term;

    return surface_density;
}

// Function computes the magnification of an image at the specified position
// for the given scaling.
//
// Input:
//      double image_position : dimensionless position of the image considered
//      double scaling_constant : characteristic scale of the profile
//
// Ouput:
//      std::complex<double> magnification : complex value of the magnification
//                                           of the specified image
std::complex<double> Magnification(double image_position,
                                   double scaling_constant) {
    double image_position_sq = image_position * image_position;
    double image_position_abs = abs(image_position);
    double image_position_abs_sq = image_position_abs * image_position_abs;

    // Compute the mass corresponding to that image
    std::complex<double> image_mass =  ImageMass(image_position_abs,
                                                 scaling_constant);

    // Compute the surface density associated with the image
    std::complex<double> surface_density = SurfaceDensity(image_position_abs,
                                                          scaling_constant);

    // Compute outer image mass term
    std::complex<double> image_mass_term = 1. - image_mass/image_position_sq;

    // Compute the bracket image mass term
    std::complex<double> image_mass_bracket_term =
        image_mass/image_position_abs_sq;

    // Compute the bracket surface density term
    std::complex<double> surface_density_bracket_term = 2. * surface_density;

    // Compute the full bracket term
    std::complex<double> bracket_term =
        1. + image_mass_bracket_term - surface_density_bracket_term;

    // Compute the magnification
    std::complex<double> magnification = 1./(image_mass_term * bracket_term);

    return magnification;
}

// Function computes the value of the morse factor for an image based on the
// dimensionless position and scaling constant. This value can take one of
// three distinct fomrs based on the sign of the surface density.
//
// Input:
//      double image_position : dimensionless position of the image considered
//      double scaling_constant : characteristic scale of the profile
//
// Output:
//      double morse_factor : value of the morse factor --- a phase change
//                            induced by the lensing
double MorseFactor(double image_position, double scaling_constant) {
    double surface_density =\
        std::real(SurfaceDensity(image_position, scaling_constant));

    double morse_factor;
    if ( surface_density == 0 ) {
        morse_factor = 0.5;
    } else if ( surface_density ) {
        morse_factor = 0.0;
    } else {
        morse_factor = 1.0;
    }

    return morse_factor;
}

// Function computes the contribution of a single image to the amplification
// factor in the geomtric optics approximation with the specified dimensionless
// frequency, source position, scaling constant, and phase.
//
// Input:
//      double image_position : dimensionless position of the image considered
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      double scaling_constant : characteristic scale of the profile
//      double phase : phase required for the minimal time delay induced by the
//                     lensing to be zero
//
// Output:
//     std::complex<double> contribution : complex value of the contribution
//                                         that the image adds to the
//                                         amplification factor
std::complex<double> ImageContribution(double image_position,
                                       double dimensionless_frequency,
                                       double source_position,
                                       double scaling_constant,
                                       double phase) {
    // i Operator
    using std::literals::complex_literals::operator""i;

    // Compute the time delay term
    std::complex<double> time_delay = TimeDelay(image_position,
                                                source_position,
                                                scaling_constant,
                                                phase);
    std::complex<double> time_delay_term =
        1i * dimensionless_frequency * time_delay;

    // Compute the Magnification term
    std::complex<double> magnification = Magnification(image_position,
                                                       scaling_constant);
    std::complex<double> magnification_term = sqrt(abs(magnification));

    // Compute the Morse Factor
    double morse_factor = MorseFactor(image_position,
                                      scaling_constant);
    std::complex<double> morse_term = 1i * M_PI * morse_factor;

    // Construct the final contribution value
    std::complex<double> contribution = exp(time_delay_term - morse_term)
                                        * magnification_term;
    
    return contribution;
}

// Function computes the amplification factor for an isolated axially symmetric
// Navarro, Frenk, White (NFW) at a given scaling using the geometric optics
// approximation for given values of dimensionless frequency and source
// position. Optionally one may give the image positions and phase if these
// have been precomputed to optimise for speed.
//
// Input:
//      double dimensionless_frequency : dimensionless form of the frequency
//                                       being amplified
//      double source_position : dimensionless displacement from optical axis
//      double scaling_constant : characteristic scale of the profile
//
// Optional Input:
//      double * image_positions : optional array of image positions
//      double phase : value of the phase required for the minimal time delay
//                     induced by the lensing to be zero
//      int number_of_images : length of the image_positions array
//
// Output:
//      std::complex<double> amplification_factor : complex value of the
//                                                  amplification factor
std::complex<double> AmplificationFactorGeometric(
    double dimensionless_frequency,
    double source_position,
    double scaling_constant,
    double * image_array = NULL,
    double phase = 0.0,
    int number_of_images = 0) {
    // i Operator
    using std::literals::complex_literals::operator""i;

    std::vector<double> image_positions;
    // Compute the Image Positions if needed
    if ( image_array != NULL && number_of_images != 0 ) {
        std::vector<double> image_position(image_array,
                                           image_array + number_of_images);
        image_positions = image_position;
    } else {
        image_positions = ImagePositions(source_position, scaling_constant);
        number_of_images = image_positions.size();
    }

    // Compute the phase for zero minimal time delay if needed
    if ( phase == 0.0 ) {
        phase = std::real(Phase(source_position, scaling_constant));
    }

    std::complex<double> amplification_factor = 0;
    // Compute the image contribution for each image
    for (int i=0; i < number_of_images; i++) {
        amplification_factor += ImageContribution(image_positions[i],
                                                  dimensionless_frequency,
                                                  source_position,
                                                  scaling_constant,
                                                  phase);
    }

    return amplification_factor;
}
