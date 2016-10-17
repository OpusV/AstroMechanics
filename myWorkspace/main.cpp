#include "orbitalMechanics.hpp"

int main(){
  
  std::cout << " --- stumpffC --- " << std::endl;
  std::cout << "stumpffC(100) " << mant::itd::stumpffC(100) << std::endl;
  std::cout << "stumpffC(-100) " << mant::itd::stumpffC(-100) << std::endl;
  std::cout << "stumpffC(0) " << mant::itd::stumpffC(0) << std::endl << std::endl;

  std::cout << " --- stumpffS --- " << std::endl;  
  std::cout << "stumpffS(100) " << mant::itd::stumpffS(100) << std::endl;
  std::cout << "stumpffS(-100) " << mant::itd::stumpffS(-100) << std::endl;
  std::cout << "stumpffS(0) " << mant::itd::stumpffS(0) << std::endl << std::endl;
  
  std::cout << " --- Lambert --- " << std::endl;
  std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> velocities = mant::itd::lambert(
    {5000.0, 10000.0, 2100.0}, 
    {-14600.0, 2500.0, 7000.0}, 
    3600.0, 
    true); 
  std::cout << "v1 " << std::endl << velocities.first << std::endl;
  std::cout << "v2 " << std::endl << velocities.second << std::endl << std::endl;
  
  std::cout << " --- Planet --- " << std::endl;
  std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> planetPosAndVel = mant::itd::positionAndVelocityOnOrbit({1.00000261, 0.01671123, -0.00001531, 0, 102.93768193, 100.46457166}, {0.00000562, -0.00004392, -0.01294668, 0, 0.32327364, 35999.37244981}, {2003,8,27,12,0,0});
  std::cout << "r " << std::endl << planetPosAndVel.first << std::endl;
  std::cout << "v " << std::endl << planetPosAndVel.second << std::endl;
  
  return 0;
}
