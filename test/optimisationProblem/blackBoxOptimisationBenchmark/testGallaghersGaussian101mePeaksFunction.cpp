// Catch
#include <catch.hpp>
#include "../../catchHelper.hpp"

class TestGallaghersGaussian101mePeaksFunction : public mant::bbob::GallaghersGaussian101mePeaksFunction {
 public:
  using mant::bbob::GallaghersGaussian101mePeaksFunction::GallaghersGaussian101mePeaksFunction;

  // Increases the visibility of internal parameters, to make them accessible.
  using mant::bbob::GallaghersGaussian101mePeaksFunction::localParameterConditionings_;
  using mant::bbob::GallaghersGaussian101mePeaksFunction::localParameterTranslations_;
  using mant::bbob::GallaghersGaussian101mePeaksFunction::rotationQ_;
};

SCENARIO("bbob::GallaghersGaussian101mePeaksFunction.getObjectiveFunctions", "[bbob::GallaghersGaussian101mePeaksFunction][bbob::GallaghersGaussian101mePeaksFunction.getObjectiveFunctions]") {
  GIVEN("A parameter") {
    THEN("Return its objective value") {
      TestGallaghersGaussian101mePeaksFunction optimisationProblem(2);

      optimisationProblem.localParameterConditionings_ = arma::mat::fixed<2, 101>({{0.177800000000000, 0.389900000000000, 0.0593000000000000, 0.0731000000000000, 0.100000000000000, 0.0870000000000000, 0.965700000000000, 0.0327000000000000, 0.480600000000000, 0.613600000000000, 0.141700000000000, 0.115000000000000, 0.256500000000000, 0.0553000000000000, 0.839900000000000, 0.239200000000000, 0.0364000000000000, 1, 0.247700000000000, 0.316200000000000, 0.168800000000000, 0.0901000000000000, 0.0756000000000000, 0.497700000000000, 0.127700000000000, 0.705500000000000, 0.515400000000000, 0.146800000000000, 0.284800000000000, 0.811100000000000, 0.0376000000000000, 0.783300000000000, 0.376500000000000, 0.0658000000000000, 0.681300000000000, 0.869700000000000, 0.294900000000000, 0.464200000000000, 0.119100000000000, 0.592600000000000, 0.136900000000000, 0.572200000000000, 0.0635000000000000, 0.635400000000000, 0.0933000000000000, 0.0966000000000000, 0.403700000000000, 0.418000000000000, 0.0705000000000000, 0.351100000000000, 0.432900000000000, 0.0351000000000000, 0.194000000000000, 0.265600000000000, 0.200900000000000, 0.0481000000000000, 0.0811000000000000, 0.448200000000000, 0.0433000000000000, 0.0498000000000000, 0.0681000000000000, 0.0418000000000000, 0.363600000000000, 0.327500000000000, 0.132200000000000, 0.223100000000000, 0.111000000000000, 0.275000000000000, 0.0464000000000000, 0.103600000000000, 0.900600000000000, 0.932600000000000, 0.174800000000000, 0.163000000000000, 0.0614000000000000, 0.187400000000000, 0.181000000000000, 0.0339000000000000, 0.730500000000000, 0.0448000000000000, 0.552600000000000, 0.0390000000000000, 0.231000000000000, 0.215400000000000, 0.157400000000000, 0.657900000000000, 0.0572000000000000, 0.0840000000000000, 0.0316000000000000, 0.0534000000000000, 0.0404000000000000, 0.339100000000000, 0.305400000000000, 0.0515000000000000, 0.107200000000000, 0.533700000000000, 0.152000000000000, 0.756500000000000, 0.0783000000000000, 0.123300000000000, 0.208100000000000}, {5.62340000000000, 2.56500000000000, 16.8761000000000, 13.6887000000000, 10.0, 11.4976000000000, 1.03550000000000, 30.5386000000000, 2.08060000000000, 1.62980000000000, 7.05480000000000, 8.69750000000000, 3.89860000000000, 18.0957000000000, 1.19060000000000, 4.18030000000000, 27.5039000000000, 1.0, 4.03700000000000, 3.16230000000000, 5.92550000000000, 11.1034000000000, 13.2194000000000, 2.00920000000000, 7.83320000000000, 1.41750000000000, 1.94030000000000, 6.81290000000000, 3.51120000000000, 1.23280000000000, 26.5609000000000, 1.27660000000000, 2.65610000000000, 15.1991000000000, 1.46780000000000, 1.14980000000000, 3.39080000000000, 2.15440000000000, 8.39930000000000, 1.68760000000000, 7.30530000000000, 1.74750000000000, 15.7387000000000, 1.57390000000000, 10.7227000000000, 10.3550000000000, 2.47710000000000, 2.39210000000000, 14.1747000000000, 2.84800000000000, 2.31010000000000, 28.4804000000000, 5.15370000000000, 3.76490000000000, 4.97700000000000, 20.8057000000000, 12.3285000000000, 2.23090000000000, 23.1013000000000, 20.0923000000000, 14.6780000000000, 23.9215000000000, 2.75040000000000, 3.05390000000000, 7.56460000000000, 4.48240000000000, 9.00630000000000, 3.63590000000000, 21.5443000000000, 9.65710000000000, 1.11030000000000, 1.07230000000000, 5.72240000000000, 6.13590000000000, 16.2975000000000, 5.33670000000000, 5.52620000000000, 29.4915000000000, 1.36890000000000, 22.3092000000000, 1.80960000000000, 25.6502000000000, 4.32880000000000, 4.64160000000000, 6.35380000000000, 1.51990000000000, 17.4753000000000, 11.9058000000000, 31.6228000000000, 18.7382000000000, 24.7708000000000, 2.94920000000000, 3.27450000000000, 19.4034000000000, 9.32600000000000, 1.87380000000000, 6.57930000000000, 1.32190000000000, 12.7662000000000, 8.11130000000000, 4.80640000000000}});
      optimisationProblem.localParameterTranslations_ = arma::mat::fixed<2, 101>({{-0.0319721223921320, 2.50613611683867, -4.70825916599147, 3.88097951076119, -1.28047655812339, 2.58250727555372, -1.54145104987895, -0.181791337107801, 0.280515878604843, -2.35738395455879, 2.26051913326145, -1.02668228157445, 1.11915203384650, 0.436659162130804, 3.40424291418613, 0.913156864010443, -1.77520187952270, 1.13488645825911, -1.08565393900114, -3.24270358859398, 2.77634326321964, -1.90113713288593, -0.481118203383186, 1.32186615230678, -4.55161493529191, 4.58647729252639, -0.113437704866684, -0.595087170770928, -1.92246892540874, -4.21803855352663, -3.77940689262812, 4.20965980223997, 2.13212318749664, -1.45071841362013, 0.963952429551355, 0.950441601838788, 0.00261979457913200, 0.129612850057174, 0.631538478137595, -0.367820698638370, -3.30261332781088, -3.01548733531064, -2.60616864280250, 0.628121124995978, -1.28365619996134, 1.25380192473267, 2.45760957272870, 3.39801939433981, -2.01923237076314, -4.82037742686823, 2.27745695178366, 1.89390394113856, 4.48779442332363, -1.98614896274293, 3.67012442456044, -2.82399720138893, 1.76652864507469, 4.26964694466519, -4.06255037812302, 0.239683218267808, 4.91506865539406, -1.76032091633217, -3.45227067556430, 2.11194049160239, -3.62817775307307, 4.25564415215959, 0.740557236664515, 0.584070133578195, -4.79518006701106, -3.54669023527563, 1.39962774174088, 0.0789440138130500, -0.532382925716046, -2.49652595499367, -2.31044045594030, -2.86985272208959, 1.41520478478538, -2.28126901593994, -2.25486047583975, 0.524013812116095, 2.55496554572812, 2.09457529412769, -0.787983837218738, 4.32266094814615, 4.76749634134431, 0.959628946851493, 0.986855957443772, 0.0371072745432470, -4.94413911824120, 0.983074860711259, 4.12552330078494, -1.76648582004984, 0.652731752579674, 2.47811352581189, -0.275333680637613, 0.677353738830577, 0.128767009538643, 1.78437111440123, 0.563872818154076, -1.72083169780272, 2.70823261413325}, {1.03733202588894, -0.0987879726864630, -4.34828945030069, -1.40182663359173, -3.75227935192573, -3.39204516236489, -2.42613681263010, 1.43687296049424, 4.31267361643233, 1.13508092488290, -2.92614086026188, 3.22913252674121, 0.0488099518192870, -0.0778442789910110, -0.723724409093243, 4.45938240539755, -3.74243000824079, 1.34067394834890, -3.14300914883592, 3.13115330290321, 1.06137205880768, 1.14932678585987, 1.03822174106562, 2.48480990819495, -0.464912095536891, -0.956787240891419, -0.813786546969722, 4.24787703263265, 0.577592379087747, 4.50035793371256, -1.23658206416477, 4.38856155691366, -0.131722731014971, 4.23626806569128, -2.11803310147908, 0.609067110613308, -0.384186889307795, -0.303887109458705, 1.61855460932072, -1.28531268952939, 4.21563603148121, -2.25081756854807, 2.84769000869571, 4.54178586508681, -4.45500849699816, -0.343789454068130, 4.01205953710770, -1.67673598495335, -2.73278123902296, 3.97078664386235, -2.22947162170183, -3.96330587391019, -2.05537756897542, 0.0170392853195720, 4.76964169257980, 4.96679502665276, -2.12056747718696, -4.69719864194744, -3.34925325758765, 0.571472759773794, -1.10765862445504, -3.46488324844539, -1.41912991357218, -2.75917649122591, 2.46204150592017, 4.94974777166868, 1.41037333194883, 3.99436842234777, -1.41706465091178, -0.675812093622085, -4.98230216563068, -1.15030786656546, -3.74170725289078, -0.0151850916878400, -0.203800888251312, 1.35372921866260, 4.97504833390335, -4.44969826505355, 4.38239202573936, 0.478619589206120, 0.984339377492162, -2.80216766450092, -2.15706713795644, 2.73734422351995, -4.13457958224210, 1.54392709097891, -2.82700569335132, 4.85801662620870, -0.211973219634888, -0.146443043799087, -2.71969075310450, -4.47956003798715, 3.49978393592325, -0.175910126050997, -2.23398116117207, 2.98904052635774, 4.95496243993079, -0.706907264981908, -2.57449244413552, 3.37269014763146, 4.37394609825683}});
      optimisationProblem.rotationQ_ = mant::rotationMatrix2d(0.1);

      CHECK(optimisationProblem.getObjectiveFunctions().at(0).first({1.0, -2.0}) == Approx(14.9835671221));
    }
  }

  THEN("Return the objective function name") {
    mant::bbob::GallaghersGaussian101mePeaksFunction optimisationProblem(2);

    CHECK(optimisationProblem.getObjectiveFunctions().size() == 1);
    CHECK(optimisationProblem.getObjectiveFunctions().at(0).second == "BBOB Gallagher's Gaussian 101-me Peaks Function (f21)");
  }
}
