#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>   // For std::inner_product or a similar concept if needed for dot product
#include <stdexcept> // For runtime_error
#include <algorithm> // For std::min and std::max

// Define a simple 3D point structure (can be replaced with your existing Point/Vector class)
struct Point3D {
    double x, y, z;

    Point3D operator-(const Point3D& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    Point3D operator+(const Point3D& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    Point3D operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    Point3D operator/(double scalar) const {
        if (scalar == 0) throw std::runtime_error("Division by zero in Point3D scalar division.");
        return {x / scalar, y / scalar, z / scalar};
    }

    double norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    Point3D normalize() const {
        double n = norm();
        if (n == 0) throw std::runtime_error("Cannot normalize a zero vector.");
        return {x / n, y / n, z / n};
    }

    static double dot(const Point3D& v1, const Point3D& v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }
};

class KappaTrajectorySolver {
public:
    // Constructor
    KappaTrajectorySolver() = default;

    /**
     * @brief Calculates the optimal kappa* value.
     *
     * @param w_prev Previous waypoint.
     * @param w_i Current waypoint (turn point).
     * @param w_next Next waypoint.
     * @param v_hat Preset speed.
     * @param c Turning rate.
     * @param eps Target precision for kappa*.
     * @return The optimal kappa* value.
     */
    double findKappaStar(const Point3D& w_prev, const Point3D& w_i, const Point3D& w_next,
                           double v_hat, double c, double eps = 1e-6) {
        if (c == 0) {
            throw std::invalid_argument("Turning rate 'c' cannot be zero.");
        }
        R_ = v_hat / c; //

        Point3D v1 = w_i - w_prev;
        Point3D v2 = w_next - w_i;

        if (v1.norm() == 0 || v2.norm() == 0) {
             // Or handle as a straight line, no turn, kappa doesn't apply or is degenerate.
            // For now, let's assume this means kappa is irrelevant, perhaps return 0 or 1
            // depending on desired behavior for straight segments.
            // For the context of the problem trying to find a specific kappa for a turn,
            // this is an invalid scenario for the Lambda function as beta would be undefined.
            // std::cerr << "Warning: Zero length segment detected. Turn angle calculation is undefined." << std::endl;
            // Based on the problem, Lambda(kappa) = 0 is the goal.
            // If there's no turn, L(kappa) = L_orig trivially.
            // However, the Lambda function itself is defined for a turn.
            // Let's assume for a straight line, any kappa is "valid" in terms of length matching,
            // but the typical use case is for turns. If it's perfectly straight (beta=pi),
            // the geometry changes.
            // The problem implies beta is an angle of a turn. If v1 or v2 is zero, beta is ill-defined.
            // If t_hat_i_ip1 and t_hat_im1_i are pointing in the same direction or opposite,
            // beta would be pi or 0.
            // Let's re-evaluate beta calculation to be robust.
            throw std::invalid_argument("Waypoints result in a zero-length segment, cannot compute turn angle.");
        }

        Point3D t_hat_im1_i = v1.normalize();
        Point3D t_hat_i_ip1 = v2.normalize();

        double dot_product = Point3D::dot(t_hat_i_ip1, t_hat_im1_i);
        // Clamp dot_product to avoid domain errors with acos due to potential floating point inaccuracies
        dot_product = std::max(-1.0, std::min(1.0, dot_product));
        beta_ = std::acos(dot_product); //

        // If beta is very close to 0 (straight segment going forward) or M_PI (straight segment reversing),
        // the turn geometry assumed by Lambda(kappa) might be degenerate.
        // The document implies beta is a turn angle, so beta != 0 and beta != PI.
        // Specifically, cot(beta/2) is used, which is problematic if beta = 0 or beta = 2*PI.
        // (pi-beta)/2 implies beta should be < pi.
        // cos(beta/2) is fine. sin(beta/2) in Xi implies beta/2 not equal to k*pi.
        if (beta_ < 1e-9 || std::abs(beta_ - M_PI) < 1e-9) {
            // For a straight line (beta = PI for forward, beta = 0 for reverse in terms of t_hat dot),
            // L(kappa) should intrinsically equal L_orig.
            // The Lambda function aims to correct deviations. If there's no turn, Lambda should be 0.
            // Let's check the components of Lambda:
            // (pi-beta)/2 -> if beta=pi, this is 0. if beta=0, this is pi/2
            // cos(beta/2) -> if beta=pi, 0. if beta=0, 1
            // cot(beta/2) -> if beta=pi, 0. if beta=0, undefined.
            // This indicates the formula is for beta in (0, pi).
            // std::cerr << "Warning: Path is nearly straight (beta approx 0 or PI). Kappa may not be well-defined by this Lambda." << std::endl;
            // If it's straight, L_kappa = L_orig for any kappa by definition as there's no turn to modify.
            // Thus Lambda(kappa) = 0. We can return any valid kappa, e.g., 0 or 1.
            // However, the premise of the bisection relies on Lambda(0) > 0 and Lambda(1) < 0 for a turn.
            // If it's perfectly straight (beta = pi), then sin(beta/2)=0, cos(beta/2)=0.
            // Xi = (1+k)/2.
            // If beta_ is very small (sharp turn, almost 180 deg U-turn):
            //  (pi-beta)/2 approaches pi/2.
            //  sin(beta/2) approaches 0. Xi approaches (1+kappa)/2.
            //  cos(beta/2) approaches 1.
            //  cot(beta/2) becomes very large. This is the critical term for small beta.
            // If beta_ is very close to PI (shallow turn, almost straight):
            //  (pi-beta)/2 approaches 0.
            //  sin(beta/2) approaches 1. Xi approaches (1+kappa + 1-kappa)/2 = 1.
            //  cos(beta/2) approaches 0.
            //  cot(beta/2) approaches 0.
            //  If Xi = 1, arccos(Xi)=0, sqrt(1/Xi^2 - 1) term is 0/0 or problematic.
            //  The document mentions numerical robustness for Xi approx 1.
            if (beta_ < 1e-6) { // Very sharp U-turn, close to 0 degrees between vectors
                 // std::cerr << "Warning: Beta is very close to 0 (sharp U-turn). Lambda function might be unstable." << std::endl;
                 // This scenario means the vehicle is nearly reversing its path.
                 // The geometry described might not be typical for such a maneuver.
                 // For now, proceed, but be aware of potential issues.
            }
            if (std::abs(beta_ - M_PI) < 1e-6) { // Almost straight
                // In this case, Lambda(kappa) should ideally be 0 for all kappa.
                // The problem states Lambda(0) > 0 and Lambda(1) < 0. This holds for turns.
                // For a straight line, the turn radius is infinite, or the concept of kappa-turn doesn't apply.
                // L(kappa) = L_orig, so Lambda(kappa) = 0. The bisection wouldn't find a unique root if f(a) and f(b) are both 0.
                // We can return a default value, e.g. 0.5, or what the DTS system expects for straight.
                // The question is about finding kappa* such that Lambda(kappa*) = 0.
                // If Lambda is always 0, any kappa works. The algorithm's premise fails.
                // However, due to floating point, it might not be exactly 0.
                // Let's assume the waypoints form a distinct turn for the algorithm to be meaningful.
                 return 0.0; // Or 1.0, or handle as per system's needs for "no effective turn".
            }
        }


        double a = 0.0; //
        double b = 1.0; //

        double lambda_a = calculateLambda(a);
        double lambda_b = calculateLambda(b);

        // Assertion from the document: Lambda(0) > 0 and Lambda(1) < 0
        // Add a small tolerance for floating point comparisons
        double check_eps = 1e-9;
        if (lambda_a < -check_eps || lambda_b > check_eps) {
            // This can happen if beta is too close to PI (almost straight line)
            // or other edge cases where the monotonicity doesn't hold as expected
            // or the values are practically zero.
            // std::cerr << "Warning: Lambda(0) or Lambda(1) does not meet expected sign conditions." << std::endl;
            // std::cerr << "  beta = " << beta_ * 180.0 / M_PI << " degrees" << std::endl;
            // std::cerr << "  Lambda(0) = " << lambda_a << std::endl;
            // std::cerr << "  Lambda(1) = " << lambda_b << std::endl;
            // If lambda_a is already ~0, kappa=0 is the solution.
            // If lambda_b is already ~0, kappa=1 is the solution.
            if (std::abs(lambda_a) < eps) return 0.0;
            if (std::abs(lambda_b) < eps) return 1.0;

            // If signs are not as expected, bisection method cannot proceed as described.
            // This might indicate an issue with the input geometry or that the assumptions
            // of Lemma 7 are not met for this specific input.
            // One possible recovery: if both are positive, means solution might be >1 (not allowed) or error.
            // If both are negative, means solution might be <0 (not allowed) or error.
            // Or, the function isn't monotonic for these params.
            // For now, throw an error if the fundamental assumption is violated significantly.
             throw std::runtime_error("Lambda function does not satisfy endpoint conditions Lambda(0)>0 and Lambda(1)<0. "
                                      "Check waypoint geometry or algorithm assumptions for this case. "
                                      "beta = " + std::to_string(beta_ * 180.0 / M_PI) + " deg, "
                                      "L(0) = " + std::to_string(lambda_a) + ", L(1) = " + std::to_string(lambda_b));
        }
        // Ensure they are strictly different for bisection to make sense.
        // If lambda_a is positive and lambda_b is negative, proceed.
        // If lambda_a is negative (or zero) and lambda_b is positive (or zero), signs are flipped from expectation.
        // The problem states Lambda(0) > 0 and Lambda(1) < 0 strictly.

        while ((b - a) > eps) { //
            double m = (a + b) / 2.0; //
            double lambda_m = calculateLambda(m);

            if (lambda_m > 0) { //
                a = m;
            } else {
                b = m;
            }
        }
        return (a + b) / 2.0; //
    }

private:
    double R_;    // Radius: v_hat / c
    double beta_; // Turn angle in radians

    /**
     * @brief Calculates the path length difference function Lambda(kappa).
     *
     * @param kappa Parameter kappa in [0,1].
     * @return Value of Lambda(kappa).
     */
    double calculateLambda(double kappa) const {
        // Xi(kappa, beta) = ((1+kappa) + (1-kappa)*sin(beta/2)) / 2 
        double sin_beta_half = std::sin(beta_ / 2.0);
        double Xi = ((1.0 + kappa) + (1.0 - kappa) * sin_beta_half) / 2.0;

        // Numerical stability for Xi close to 0 or 1
        // Using a small delta, e.g., 1e-9, to prevent domain errors with acos or sqrt
        const double delta = 1e-12; // Adjusted delta for robustness
        Xi = std::max(delta, std::min(Xi, 1.0 - delta));
        // If Xi is exactly 1 (e.g. beta=pi and kappa=any, or beta<pi and kappa=1 leads to sin_beta_half=1),
        // then sqrt(1/Xi^2 - 1) will be sqrt(1-1)=0.
        // If Xi is exactly 0 (not possible given its formula if beta/2 in [0, pi/2] and kappa in [0,1]),
        // then 1/Xi^2 is problematic. Smallest Xi happens when kappa=1 and sin(beta/2)=0 (beta=0), so Xi=1.
        // Or kappa=-1 (not allowed) and sin(beta/2)=1 (beta=pi), then Xi = 0.
        // Given kappa in [0,1] and beta in (0,pi):
        // min sin_beta_half is close to 0 (for beta close to 0)
        // max sin_beta_half is 1 (for beta = pi, but we exclude beta=pi based on cot)
        // If beta is close to 0: sin_beta_half ~ 0. Xi ~ (1+kappa)/2. Smallest is 0.5 (kappa=0). Largest 1 (kappa=1).
        // If beta is pi (excluded, but for limit): sin_beta_half=1. Xi = ((1+k)+(1-k)*1)/2 = (1+k+1-k)/2 = 1.
        // So Xi should be in [~0.5, 1]. Clamping to 1-delta is most important.

        double term_acos_Xi = std::acos(Xi);
        double term_sqrt_Xi;
        if (std::abs(Xi) < 1.0 - 1e-9) { // Avoid division by zero if Xi is 1.0
             term_sqrt_Xi = Xi * std::sqrt(1.0 / (Xi * Xi) - 1.0);
        } else if (Xi >= 1.0 - 1e-9) { // If Xi is very close to 1 (or 1)
            term_sqrt_Xi = 0.0; // Limit of Xi * sqrt(1/Xi^2 - 1) as Xi -> 1 is 0
            term_acos_Xi = 0.0; // acos(1) = 0
        }
        else { // Xi is very close to -1 (should not happen with kappa in [0,1] and beta in [0,PI])
            term_sqrt_Xi = 0.0; // sqrt would be complex
            term_acos_Xi = M_PI; // acos(-1) = PI
        }


        double cos_beta_half = std::cos(beta_ / 2.0);
        double cot_beta_half;
        if (std::abs(std::sin(beta_ / 2.0)) < 1e-9) { // beta is close to 0 or 2*PI
            // This case should ideally be handled by the beta check in findKappaStar
            // If beta_ is very small, cot_beta_half becomes huge.
            // Based on prior checks, beta_ should not be 0.
            throw std::runtime_error("cot(beta/2) undefined: beta is too close to 0 or 2*PI.");
        } else {
            cot_beta_half = cos_beta_half / sin_beta_half;
        }

        /* Original formula:
           Lambda(kappa) = 2R * ( (pi-beta)/2
                                 + 2*acos(Xi)
                                 - Xi*sqrt(1/Xi^2 - 1)
                                 - (1-kappa)*cos(beta/2)
                                 - kappa*cot(beta/2) )
        */
        double lambda_val = (M_PI - beta_) / 2.0
                          + 2.0 * term_acos_Xi
                          - term_sqrt_Xi
                          - (1.0 - kappa) * cos_beta_half
                          - kappa * cot_beta_half;

        return 2.0 * R_ * lambda_val;
    }

#ifndef M_PI
    static constexpr double M_PI = 3.14159265358979323846;
#endif
};

// --- Example Usage ---
// To compile (e.g., with g++):
// g++ -std=c++17 your_file_name.cpp -o kappa_solver
// ./kappa_solver
int main() {
    KappaTrajectorySolver solver;

    // Example Test Case (you should replace these with actual values from your simulation)
    Point3D w_prev = {0, 0, 0};
    Point3D w_i    = {5, 0, 0};
    Point3D w_next = {5, 5, 0}; // 90 degree turn

    double v_hat = 10.0; // m/s
    double c     = 1.0;  // rad/s (turn rate)
    double eps   = 1e-7;

    std::cout << "Test Case 1: 90-degree turn" << std::endl;
    try {
        double kappa_star = solver.findKappaStar(w_prev, w_i, w_next, v_hat, c, eps);
        std::cout << "  Optimal kappa*: " << kappa_star << std::endl;

        // Verify Lambda(kappa_star) is close to 0
        // Need to re-setup solver's internal R_ and beta_ for an external check,
        // or make calculateLambda public if direct testing is needed.
        // For simplicity, we'll trust the internal calculation for now.
        // To test lambda(kappa_star) ~ 0, one would need to call calculateLambda directly.
        // This implies calculateLambda might need to be public or a helper test function.
    } catch (const std::exception& e) {
        std::cerr << "  Error: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // Test Case 2: Sharper turn
    Point3D w_next_sharp = {6, 1, 0}; // Approx 45 deg turn angle with (5,0)-(0,0) and (6,1)-(5,0)
                                      // Vector 1: (5,0,0). Vector 2: (1,1,0)
                                      // Angle between (5,0) and (1,1) from origin: (1,1) is at 45 deg. Vector from (5,0) to (6,1) is (1,1)
                                      // So turn from x-axis to vector (1,1) means beta is 135 deg for the path vectors.
                                      // v_prev_i = (5,0,0), v_i_next = (1,1,0)
                                      // acos( (5,0,0) dot (1,1,0) / (|5|*|sqrt(2)|) ) = acos(5 / (5*sqrt(2))) = acos(1/sqrt(2)) = 45 deg. This is beta.
                                      // The document defines beta = acos( t_hat_i_ip1^T * t_hat_im1_i )
                                      // t_im1_i = normalize(w_i - w_prev) = normalize( (5,0,0) - (0,0,0) ) = (1,0,0)
                                      // t_i_ip1 = normalize(w_next - w_i) = normalize( (6,1,0) - (5,0,0) ) = normalize( (1,1,0) ) = (1/sqrt(2), 1/sqrt(2), 0)
                                      // dot_product = 1 * 1/sqrt(2) + 0 * 1/sqrt(2) = 1/sqrt(2)
                                      // beta = acos(1/sqrt(2)) = pi/4 = 45 degrees.

    w_prev = {0,0,0};
    w_i    = {5,0,0};
    w_next = {5+cos(M_PI/4.0), 5+sin(M_PI/4.0), 0}; // Should be w_i + vector, e.g. {5+1, 1, 0} for a clear next step.
                                               // Let's define w_next to make a clear angle.
                                               // w_next such that (w_next - w_i) makes 45 deg with (w_i - w_prev)
                                               // w_i - w_prev = (5,0,0)
                                               // Let w_next - w_i = (1,1,0). So w_next = (6,1,0)
    w_next = {6,1,0};


    std::cout << "Test Case 2: 45-degree turn (beta = 45 deg)" << std::endl;
    try {
        double kappa_star_sharp = solver.findKappaStar(w_prev, w_i, w_next, v_hat, c, eps);
        std::cout << "  Optimal kappa*: " << kappa_star_sharp << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  Error: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // Test Case 3: Shallow turn (beta close to PI)
    w_prev = {0,0,0};
    w_i    = {5,0,0};
    w_next = {10, 0.1, 0}; // Very shallow turn
                           // t_im1_i = (1,0,0)
                           // t_i_ip1 = normalize( (5, 0.1, 0) ). Dot product very close to 1. Beta close to 0 for these vectors.
                           // The definition of beta is acos(t_hat_i_ip1^T * t_hat_im1_i)
                           // Here, t_hat_im1_i is direction of (w_i - w_prev)
                           // t_hat_i_ip1 is direction of (w_next - w_i)
                           // If w_prev=(0,0,0), w_i=(5,0,0), w_next=(10,0.1,0):
                           // t_hat_im1_i = (1,0,0)
                           // vec_i_ip1 = (5, 0.1, 0). norm = sqrt(25.01) approx 5.001
                           // t_hat_i_ip1 approx (5/5.001, 0.1/5.001, 0) approx (0.9998, 0.0199, 0)
                           // dot_product approx 0.9998. acos(0.9998) approx 0.02 radians (1.14 degrees). This is beta.
                           // This is a very SHARP turn in terms of the algorithm's beta (angle between successive segment vectors).
                           // The problem says: β = cos⁻¹(t̂ᵢ,ᵢ₊₁ᵀ t̂ᵢ₋₁,ᵢ)
                           // For my current setup: w_prev -> w_i is along +x.
                           // w_i -> w_next is a vector (5, 0.1, 0).
                           // The angle between (1,0,0) and (5,0.1,0) is small.
                           // This means beta is small.
                           // A "shallow turn" in vehicle path means the angle of deviation from straight is small.
                           // If beta is the internal angle of the corner (like in a polygon), then a shallow turn has beta close to PI.
                           // The document refers to (pi-beta)/2. If beta is the turn angle (0 for straight, pi for U-turn),
                           // then (pi-beta)/2 means beta is the angle of deviation from straight, not the interior angle.
                           // The pseudocode uses: β = arccos(dot(normalize(w_next - w_i), normalize(w_i - w_prev)))
                           // This is the angle between the two path segments.
                           // If w_i-w_prev = (1,0) and w_next-w_i = (1,0.01), beta is small. This is a sharp deviation *from the previous segment's direction*.
                           // This is what the code calculates. A "shallow turn" in terms of path curvature would mean beta is close to PI (180 degrees),
                           // meaning the change in direction is small.
                           // Let's re-evaluate the "shallow turn" example:
                           // w_prev = (0,0,0), w_i = (5,0,0), w_next = (4,1,0) -> results in beta > 90 deg
                           // This means t_im1_i is (1,0,0) and t_i_ip1 is normalize(-1,1,0). dot is -1/sqrt(2). beta = 135 deg. This is a "gentle" turn.
    w_next = {4, 1, 0}; // beta = 135 degrees (gentle turn)
    std::cout << "Test Case 3: Gentle turn (beta = 135 deg)" << std::endl;
    try {
        double kappa_star_shallow = solver.findKappaStar(w_prev, w_i, w_next, v_hat, c, eps);
        std::cout << "  Optimal kappa*: " << kappa_star_shallow << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  Error: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;


     // Test Case 4: Almost straight line (beta close to PI), where Lambda(0) and Lambda(1) might be problematic
    w_prev = {0,0,0};
    w_i    = {5,0,0};
    w_next = {10,0.00001,0}; // Almost straight - (w_next - w_i) is almost same direction as (w_i - w_prev)
                             // t_im1_i = (1,0,0)
                             // t_i_ip1 = normalize(5, 0.00001, 0) is very close to (1,0,0).
                             // dot product is very close to 1. beta is very close to 0.
                             // This is a *very sharp* turn according to the beta definition (vehicle almost continues straight).
                             // No, this is confusing. Beta = arccos(dot(v2,v1)).
                             // If v1 = (1,0) and v2 = (1,epsilon), then dot is close to 1, beta is close to 0.
                             // This is a very *slight* turn. A shallow turn where the change in heading is small.
                             // The document's (pi-beta)/2 implies beta is the *interior* angle of the turn corner if we view segments.
                             // Let's re-read "β = cos⁻¹(t̂ᵢ,ᵢ₊₁ᵀ t̂ᵢ₋₁,ᵢ)"
                             // t_hat_im1_i points from w_im1 to w_i.
                             // t_hat_i_ip1 points from w_i to w_ip1.
                             // If straight: t_hat_i_ip1 = t_hat_im1_i. Dot product is 1. Beta = 0.
                             // If 180 deg U-turn: t_hat_i_ip1 = -t_hat_im1_i. Dot product is -1. Beta = PI.
                             // So beta IS the turn angle: 0 for straight, PI for full reversal.
                             // My previous interpretation was mixed up.
                             // (pi-beta)/2:
                             //  - If straight (beta=0), then pi/2.
                             //  - If 90 deg turn (beta=pi/2), then pi/4.
                             //  - If U-turn (beta=pi), then 0.
                             // cot(beta/2):
                             //  - If straight (beta=0), then cot(0) is undefined. This is an issue.
                             //  - If 90 deg turn (beta=pi/2), then cot(pi/4) = 1.
                             //  - If U-turn (beta=pi), then cot(pi/2) = 0.
                             // This means the formula has a singularity or problem at beta=0 (straight).
                             // The document states Lambda(0)>0, Lambda(1)<0 (Lemma 7).
                             // This implies beta is not 0. "Thus, there exists ... unique k* in (0,1)" (Theorem 8).
                             // The problem is defined for a *turn*. If it's straight, L(kappa)=L_orig, so Lambda=0 for all kappa.
                             // The algorithm should handle this. My `beta_ < 1e-6` check returning 0.0 is one way.
    std::cout << "Test Case 4: Almost straight (beta close to 0)" << std::endl;
    try {
        // For beta close to 0, the initial check `if (beta_ < 1e-6)` should catch it.
        // If it's truly straight, my code might return 0.0 from that early exit.
        // The `findKappaStar` will throw if `cot(beta/2)` is problematic due to `beta_` near 0.
        // My current `beta_ < 1e-6` check returns 0.0 for kappa*.
        // Let's make it slightly non-zero beta
        w_next = {10, 0.001, 0}; // beta will be small but non-zero.
        double kappa_star_straightish = solver.findKappaStar(w_prev, w_i, w_next, v_hat, c, eps);
        std::cout << "  Optimal kappa*: " << kappa_star_straightish << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  Error for 'almost straight': " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // Test Case 5: U-turn (beta close to PI)
    w_prev = {0,0,0};
    w_i    = {5,0,0};
    w_next = {0,0.001,0}; // Almost a U-turn relative to w_i. (w_next - w_i) is (-5, 0.001, 0)
                          // t_im1_i = (1,0,0)
                          // t_i_ip1 = normalize(-5, 0.001, 0), approx (-1, small_val, 0)
                          // Dot product close to -1. Beta close to PI.
    std::cout << "Test Case 5: Almost U-turn (beta close to PI)" << std::endl;
    try {
        // My `std::abs(beta_ - M_PI) < 1e-6` check returns 0.0 for kappa*.
        // Let's make beta slightly less than PI
        w_next = {0, 0.1, 0}; // beta will be large, close to PI.
        double kappa_star_uturn = solver.findKappaStar(w_prev, w_i, w_next, v_hat, c, eps);
        std::cout << "  Optimal kappa*: " << kappa_star_uturn << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  Error for 'almost U-turn': " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;


    return 0;
}