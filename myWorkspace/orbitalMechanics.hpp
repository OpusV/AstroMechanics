#pragma once

// C++ standard library
#include <cmath>
#include <iostream>
#include <functional>

// Armadillo
#include <armadillo>

// Mantella
#include <mantella>

namespace mant{
  namespace itd{
    
    double stumpffC(
        const double parameter){
      if(parameter > 0){
        return (1.0 - std::cos(std::sqrt(parameter))) / parameter;
      }else if (parameter < 0){
        return (std::cosh(std::sqrt(-parameter)) - 1.0) / -parameter;
      }else{
        return 0.5;
      }
    }
    
    double stumpffS(
        const double parameter){
      
      if(parameter > 0){
        return (std::sqrt(parameter) - std::sin(std::sqrt(parameter))) / std::pow(std::sqrt(parameter), 3.0);
      }else if (parameter < 0){
        return (std::sinh(std::sqrt(-parameter)) - std::sqrt(-parameter)) / std::pow(std::sqrt(-parameter), 3.0);
      }else{       
        return 1.0 / 6.0;  
      }
    }
    
    std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> lambert(
        const arma::Col<double>::fixed<3> &startPosition,
        const arma::Col<double>::fixed<3> &endPosition,
        const double flightTime,
        const bool isPrograde){
      
      double startPositionNorm = arma::norm(startPosition);
      double endPositionNorm = arma::norm(endPosition);
      
      arma::Col<double> positionsCrossProduct = arma::cross(startPosition, endPosition);
      double theta = std::acos(arma::norm_dot(startPosition, endPosition));
      
      if(isPrograde){
        if(positionsCrossProduct(2) < 0){
          theta = 2.0 * arma::datum::pi - theta;
        }
      }else{
        if(positionsCrossProduct(2) >= 0){
          theta = 2.0 * arma::datum::pi - theta;
        }
      }

      double a = std::sin(theta) * std::sqrt((startPositionNorm * endPositionNorm) / (1.0 - std::cos(theta)));
      
      std::function<double(const double)> yFunction = [&](
          const double parameter){ 
        return startPositionNorm + endPositionNorm + a * ((parameter * stumpffS(parameter) - 1.0) / std::sqrt(stumpffC(parameter)));
      };
      
      double z = mant::brent(
        [&](
            const double parameter) { 

          return stumpffS(parameter) * std::pow(yFunction(parameter) / stumpffC(parameter), 3.0/2.0) + a * std::sqrt(yFunction(parameter)) - std::sqrt(3.986e+5) * flightTime; //heliocentric should be 1.32712440018e11
          
        }, -2.0, 2.0, 100); //TODO: remove magic numbers
      
      double y = yFunction(z);

      double f = 1.0 - y / startPositionNorm;
      double g = a * std::sqrt(y / 3.986e+5); //heliocentric should be 1.32712440018e11
      double gDot = 1.0 - y / endPositionNorm;
      
      return {{(1.0 / g) * (endPosition - f * startPosition)}, {(1.0 / g) * (gDot * endPosition - startPosition)}};
      
    }
    
    
    std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> positionAndVelocityOnOrbit(
        const arma::Col<double>::fixed<6>& keplerValues,
        const arma::Col<double>::fixed<6>& dKeplerValues,
        const arma::Col<double>::fixed<6>& date) {
      // parameters
      // -keplerValues as known (a, e, i, omega, w, L)
      // -date in [yyyy,(mon)mon, (d)d, (h)h, (min)min, (s)s]
      
      double jd = 367.0 * date(0) - std::trunc((7.0 * (date(0) + std::trunc((date(1) + 9.0) / 12.0))) / 4.0) + std::trunc(275.0 * date(1) / 9.0) + date(2) + 1721013.5 + (date(3) + date(4) / 60.0 + date(5) / 3600.0) / 24;
      
      double t0 = (jd - 2451545.0) / 36525.0;  
      
      arma::Col<double>::fixed<6> keplerValAtTime = keplerValues + dKeplerValues * t0;
      keplerValAtTime(0) = keplerValAtTime(0) * 149597871.464; //au to km
      keplerValAtTime(2) = std::fmod(keplerValAtTime(2), 360.0) * arma::datum::pi / 180.0;
      keplerValAtTime(3) = std::fmod(keplerValAtTime(3), 360.0) * arma::datum::pi / 180.0;
      keplerValAtTime(4) = std::fmod(keplerValAtTime(4), 360.0) * arma::datum::pi / 180.0;
      keplerValAtTime(5) = std::fmod(keplerValAtTime(5), 360.0) * arma::datum::pi / 180.0;   
      
      double angularMomentum = std::sqrt(1.32712440018e11 * keplerValAtTime(0) * (1.0 - std::pow(keplerValAtTime(1), 2.0))); //sun gravity!!!
      
      double argumentPerihelion = keplerValAtTime(4) - keplerValAtTime(3);
      double meanAnomaly = keplerValAtTime(5) - keplerValAtTime(4); 
      
      double eccentricAnomaly0;
      if(meanAnomaly < arma::datum::pi) {       
        eccentricAnomaly0 =  meanAnomaly + keplerValAtTime(1) / 2.0;       
      } else {        
        eccentricAnomaly0 =  meanAnomaly - keplerValAtTime(1) / 2.0;        
      }
      
      double eccentricAnomaly = mant::brent(
        [&](
            const double parameter) {

          return parameter - keplerValAtTime(1) * std::sin(parameter) - meanAnomaly;
          
        }, -10.0, 10.0, 100); //TODO: remove magic numbers  
      
      double trueAnomaly = std::fmod(2.0 * std::atan(std::sqrt((1.0 + keplerValAtTime(1)) / (1.0 - keplerValAtTime(1))) * std::tan(eccentricAnomaly / 2.0)), 2.0 * arma::datum::pi);
      
      arma::Col<double>::fixed<3> helperVectorPos = {std::cos(trueAnomaly), std::sin(trueAnomaly), 0};
      arma::Col<double>::fixed<3> perifocalPosition = (std::pow(angularMomentum, 2.0) / 1.32712440018e11) * (1.0 / (1.0 + keplerValAtTime(1) * std::cos(trueAnomaly))) * helperVectorPos; //sun gravity!!!
      
      arma::Col<double>::fixed<3> helperVectorVel = {-std::sin(trueAnomaly), (keplerValAtTime(1) + std::cos(trueAnomaly)), 0};
      arma::Col<double>::fixed<3> perifocalVelocity = (1.32712440018e11 / angularMomentum) * helperVectorVel; //sun gravity!!!    
      
      arma::Mat<double>::fixed<3, 3> R3_W = {
        {std::cos(keplerValAtTime(3)),  std::sin(keplerValAtTime(3)), 0.0}, 
        {-std::sin(keplerValAtTime(3)), std::cos(keplerValAtTime(3)), 0.0}, 
        {0.0,                           0.0,                          1.0}
      };
      
      arma::Mat<double>::fixed<3, 3> R1_i = {
        {1.0, 0.0, 0.0}, 
        {0.0, std::cos(keplerValAtTime(2)), std::sin(keplerValAtTime(2))}, 
        {0.0, std::sin(keplerValAtTime(2)), std::cos(keplerValAtTime(2))}
      };
      
      arma::Mat<double>::fixed<3, 3> R3_w = {
        {std::cos(argumentPerihelion),  std::sin(argumentPerihelion), 0.0}, 
        {-std::sin(argumentPerihelion), std::cos(argumentPerihelion), 0.0}, 
        {0.0,                           0.0,                          1.0}
      };
      
      arma::Mat<double>::fixed<3, 3> Q_pX = (R3_w * R1_i * R3_W).t();
      
      return {{Q_pX * perifocalPosition},{Q_pX * perifocalVelocity}};
    }
    
  }
}
