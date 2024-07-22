#include <Rcpp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

// [[Rcpp::plugins("cpp11")]]
using namespace Rcpp;

class infected_man {
public:
  int mortality;
  std::vector<int> infection_history;
  const double weight;
  std::vector<int> five_year_infection_vec;
  int positive_on_navdx;
  int positive_on_e6;
  int positive_on_hpv;
  int negative_on_hpv;
  int screen_eligible;

  int age_in_2021;

  infected_man(const Rcpp::NumericVector& input_data,
	       const Rcpp::NumericMatrix& weights)
    : infection_history(70, 0),
      weight{weights(input_data[1] - 15, input_data[0])},
      mortality{static_cast<int>(input_data[5])},
      five_year_infection{0},
      five_year_infection_vec(70, 0),
      positive_on_navdx{0},
      positive_on_e6{0},
      positive_on_hpv{0},
      screen_eligible{1},
      negative_on_hpv{0},
      age_in_2021{static_cast<int>(input_data[1])}
  {
    double r = R::rnorm(13, 5 / 1.96);
    if (r < 1)
      r = 1;

    for (int i = 0; i < 70; i++) {
      infection_history[i] = input_data[i + 6];
      if (five_year_infection || input_data[i+6] > r) {
	five_year_infection_vec[i] = 1;
	five_year_infection = 1;
      }
    }
  };
  
private:
  int five_year_infection;
};
  
class stage_man {
public:
  int mortality_age_no_cancer;
  int mortality_age_clinical;
  int mortality_age_screen;
  std::vector<int> infection_history;
  const double weight;
  const int diagnosis_age;
  int diagnostic_stage = 0; // ajcc7 stage
  std::vector<double> stage_ages;
  std::vector<int> T;
  std::vector<int> N;
  std::vector<double> tumor_size;
  std::vector<double> node_size;
  std::vector<int> five_year_infection_vec;
  std::vector<int> ajcc8_t;
  std::vector<int> ajcc8_n;
  std::vector<int> ajcc8_stage;
  int positive_on_navdx;
  int positive_on_e6;
  int positive_on_hpv;
  int negative_on_hpv;
  int screen_eligible;
  const int causal_infection_age;
  int smk;
  int age_in_2021;
  int screen_death;
  int clinical_death;

  stage_man(const Rcpp::NumericVector& input_data,
	    const Rcpp::NumericMatrix& weights,
	    const Rcpp::NumericMatrix& mortality_data,
	    const Rcpp::NumericMatrix& stage_probs,
	    const Rcpp::NumericMatrix& progression_probs,
	    const int progression_speed,
	    const Rcpp::List& TNM)
    : infection_history(70, 0),
      weight{weights(input_data[1] - 15, input_data[0])},
      causal_infection_age{static_cast<int>(input_data[3])},
      diagnosis_age{static_cast<int>(input_data[4])},
      stage_ages(4,0),
      T(4, 0),
      N(4,0),
      tumor_size(4, 0),
      node_size(4, 0),
      ajcc8_t(4, 0),
      ajcc8_n(4, 0),
      ajcc8_stage(4, 0),
      five_year_infection{0},
      five_year_infection_vec(70, 0),
      positive_on_navdx{0},
      positive_on_e6{0},
      positive_on_hpv{0},
      negative_on_hpv{0},
      screen_eligible{1},
      smk(0),
      mortality_age_clinical(0),
      mortality_age_screen(0),
      age_in_2021{static_cast<int>(input_data[1])},
      screen_death(0),
      clinical_death(0)
  {

    double r = R::rnorm(13, 5 / 1.96);
    if (r < 1)
      r = 1;

    for (int i = 0; i < 70; i++) {
      infection_history[i] = input_data[i + 6];
      if (five_year_infection || input_data[i+6] > r) {
	five_year_infection_vec[i] = 1;
	five_year_infection = 1;
      }
    }

    if (input_data[4] < 69) {
	for (int i = input_data[4] + 1; i < 70; i++) {
		mortality_age_no_cancer = i;
	  if (R::runif(0.0, 1.0) < mortality_data(i, input_data[0])) {
		break;
	  }
	}
    } else {
      mortality_age_no_cancer = 69;
    }

    assign_diagnostic_stage(stage_probs);
    assign_TN(TNM);
    assign_prior_stages(progression_probs, progression_speed - 1);
    std::vector<int> smk_inds{4, 5, 6, 7, 12, 13, 14, 15, 20, 21, 22, 23, 28, 29, 30, 31};
    if (std::find(smk_inds.begin(), smk_inds.end(), input_data[0]) != smk_inds.end()){
      smk = 1;
    }
    
  };

  void assign_diagnostic_stage(const Rcpp::NumericMatrix& stage_probs)
  {
    double r = R::runif(0.0, 1.0);
    for (int i  = 0; i < 3; i++) {
      if (r > stage_probs(diagnosis_age, i))
	diagnostic_stage++;
    }
  };

  void assign_prior_stages(const Rcpp::NumericMatrix& progression_probs,
			   const int progression_speed)
  {
    std::vector<double> tmp_stages(4);
    tmp_stages[0] = R::rexp(1 / progression_probs(progression_speed, 1));
    
    if (diagnostic_stage > 0) {
	tmp_stages[1] = R::rexp(1 / progression_probs(progression_speed, 3));
    }
    if (diagnostic_stage > 1) {
	tmp_stages[2] = R::rexp(1 / progression_probs(progression_speed, 5));
    }
    if (diagnostic_stage > 2) {
	tmp_stages[3] = R::rexp(1 / progression_probs(progression_speed, 7));
    }
    
    tmp_stages[diagnostic_stage] = R::runif(0.0, tmp_stages[diagnostic_stage]);

    if (T[2] == 1)
      tmp_stages[1] = 0;
    
    double pre_clinical_time = 0;

    for (int i = 0; i < 4; i++)
      pre_clinical_time += tmp_stages[i];

    double time_from_infection = diagnosis_age - causal_infection_age;

    // there must be at least some time from infection to cancer
    if (time_from_infection == 0)
      time_from_infection = 1;

    // Check if cancer time exceeds infection length
    if (pre_clinical_time > time_from_infection) {
      for (int i = 0; i < 4; i++)
	tmp_stages[i] /= time_from_infection;
    }

      stage_ages[diagnostic_stage] = diagnosis_age - tmp_stages[diagnostic_stage];
      for (int i = diagnostic_stage - 1; i > -1; i--) {
	stage_ages[i] = stage_ages[i + 1] - tmp_stages[i];
      }
  };

  void assign_TN(const Rcpp::List& TNM)
  {
    int index;
    int number_options;
    int first = 1;
    for (int i = 3; i > -1; i--) {
      if (diagnostic_stage < i)
	continue;
      if (diagnostic_stage > i)
	first = 0;
      Rcpp::NumericMatrix data = TNM[i];
      number_options = data.nrow();

      while(1) {
      index = 0;
	double r = R::runif(0.0, 1.0);
	for (int j = 0; j < number_options; j++) {
	    if (r > data(j, 3)) {
	    index++;
	    } else {
	    break;
	    }
	}
	T[i] = data(index, 0);
	N[i] = data(index, 1);
	ajcc8_stage[i] = data(index, 8) - 1;
	

	double r_1 = R::runif(0.0, 1.0);
	double r_2 = R::runif(0.0, 1.0);
	if (first == 1 || tumor_size[i + 1] > 0) {
	  if (r_1 < data(index, 13)) {
	      tumor_size[i] = 0;
	    }
	  else if (tumor_size[i + 1] == 1 || r_1 < data(index, 13) + data(index, 14)) {
	    tumor_size[i] = 1;
	  } else {
	    tumor_size[i] = 2;
	  }
	}

	if (data(index, 5)) {
	  node_size[i] = 2;
	} else if (data(index, 6)) {
	  node_size[i] = 1;
	} else if (data(index, 7)) {
	  node_size[i] = 0;
	}
	if (first)
	  break;
	if ((N[i] <= N[i + 1] && T[i] <= T[i + 1] && node_size[i] <= node_size[i + 1]) ||
	    (i == 1 && T[i + 1] == 1 && node_size[i] <= node_size[i + 1])) {
	  break;
	}
      }
    }
  };
  
private:
  int five_year_infection;
};

class screen_data {
public:
  Rcpp::NumericMatrix screen_age_stage_ajcc8;
  Rcpp::NumericMatrix tx_by_age;
  std::vector<double> number_of_infection_screens;
  std::vector<double> number_of_cancer_screens;
  std::vector<double> cancer_rate_denom;
  std::vector<double> cancer_rate_num;
  Rcpp::NumericMatrix screen_age_stage_ajcc7;
  Rcpp::NumericMatrix life_results_1;
  Rcpp::NumericMatrix life_results_2;
  int inf_screens;
  int canc_screens;
  int sensitivity_2;
  int sensitivity_3;

  screen_data(int n,
	      int sensitivity_2,
	      int sensitivity_3)
    :screen_age_stage_ajcc8(70, 8),
     tx_by_age(70, 4),
     screen_age_stage_ajcc7(70, 8),
     number_of_infection_screens(8, 0),
     number_of_cancer_screens(4, 0),
     cancer_rate_denom(70, 0),
     cancer_rate_num(70, 0),
     life_results_1(n, 13),
     life_results_2(n, 13),
     inf_screens(0),
     canc_screens(0)
  {
  };


  // test 1: oral HPV16
  void screen_test_1(stage_man& staged_man,
		     int age,
		     int& need_cancer_screen)
  {
    int hpv = staged_man.infection_history[age] > 0;
    int cancer = staged_man.stage_ages[0] <= age;

    int test_result = 0;
    if (hpv == 1) {
	double r = R::runif(0.0, 1.0);
	if (r < 0.81)
	  test_result = 1;
      }

    if (test_result) {
      staged_man.positive_on_hpv = 1;
    } else {
      staged_man.positive_on_hpv = 0;
    }
  };

    void screen_test_1(infected_man& staged_man,
		       int age,
		       int& need_cancer_screen)
  {
    int hpv = staged_man.infection_history[age] > 0;

    if (hpv == 1 && (R::runif(0.0, 1.0) < 0.81)) {
      staged_man.positive_on_hpv = 1;
    } else {
      staged_man.positive_on_hpv = 0;
    }
  };

 // test 1: E6
 void screen_test_2(stage_man& staged_man,
		     int age,
		     int& need_cancer_screen)
  {
    int cancer = staged_man.stage_ages[0] <= age;
    int five_year_infection = staged_man.five_year_infection_vec[age];

    int test_result = 0;

    double r = R::runif(0.0, 1.0);
    if ((five_year_infection && r < 0.92) ||
	(!five_year_infection && cancer && r < 0.92)) {
      test_result = 1;
    }
      

    if (test_result) {
      staged_man.positive_on_e6 = 1;
    }
  };

    void screen_test_2(infected_man& staged_man,
		       int age,
		       int& need_cancer_screen)
  {
    int five_year_infection = staged_man.five_year_infection_vec[age];
    
    double r = R::runif(0.0, 1.0);
    if (five_year_infection  && r < 0.92) {
      staged_man.positive_on_e6 = 1;
    }
  };

  void screen_test_3(stage_man& staged_man,
		     int age,
		     int& need_cancer_screen)
  {
    int stage = -1;
    for (int j = 0; j < 4; j++) {
      if (staged_man.stage_ages[j] <= age && staged_man.stage_ages[j] != 0)
	stage++;
    }

    double r = R::runif(0.0, 1.0);
    int test_result = 0;

    if (staged_man.positive_on_navdx) {
      test_result = 1;
    } else if (
	       (staged_man.N[stage] == 0 && staged_man.tumor_size[stage] == 0 && r < 0.36) ||
	       (staged_man.N[stage] == 0 && staged_man.tumor_size[stage] == 1 && r < 0.36) ||
               (staged_man.N[stage] == 0 && staged_man.tumor_size[stage] == 2 && r < 0.36) ||
	       (staged_man.N[stage] >= 1 && r < 0.36) ||
		(staged_man.N[stage] == 1 && r < 0.94) ||
               (staged_man.N[stage] == 2 && r < 0.96) ||
	       (staged_man.N[stage] == 3 && r < 1)) {
      test_result = 1;
      staged_man.positive_on_navdx = 1;
    }

  };

  int get_mort_age_diff(int diagnosis_age,
			stage_man& staged_man,
			Rcpp::NumericMatrix prob_mat_screen,
			Rcpp::NumericMatrix prob_mat_clinical)
  {
    Rcpp::NumericVector probs_clinical = prob_mat_clinical(staged_man.diagnosis_age, _);
    
    double r = R::runif(0.0, 1.0);
    int years_to_clinical_death = 1;
    for (int i = 0; i < 5; i++) {
      if (r < 1 - probs_clinical[i]) {
	break;
      } else {
	years_to_clinical_death++;
      }
    }


    if (years_to_clinical_death > 5 ||
	staged_man.diagnosis_age + years_to_clinical_death >= staged_man.mortality_age_no_cancer) {
      staged_man.mortality_age_clinical = staged_man.mortality_age_no_cancer;
    } else {
      staged_man.mortality_age_clinical = staged_man.diagnosis_age + years_to_clinical_death;
      staged_man.clinical_death = 1;
    }

    if (diagnosis_age < staged_man.diagnosis_age) {
	Rcpp::NumericVector probs_screen = prob_mat_screen(diagnosis_age, _);

	double r = R::runif(0.0, 1.0);
	int years_to_screen_death = 1;
	if (years_to_clinical_death < 6 || sensitivity_3) {
	    for (int i = 0; i < 5; i++) {
		if (r < 1 - probs_screen[i]) {
		    break;
		} else {
		    years_to_screen_death++;
		}
	    }
	} else {
	  years_to_screen_death = 6;
	}

	if (years_to_screen_death > 5 ||
	    diagnosis_age + years_to_screen_death >= staged_man.mortality_age_no_cancer) {
	    staged_man.mortality_age_screen = staged_man.mortality_age_no_cancer;
	} else {
	    staged_man.mortality_age_screen = diagnosis_age + years_to_screen_death;
	    staged_man.screen_death = 1;
	}

    } else {
      staged_man.mortality_age_screen = staged_man.mortality_age_clinical;
    }

    // Man may die from other causes before clinical cancer
    if (staged_man.clinical_death == 0 && !sensitivity_3)
      staged_man.screen_death = 0;

    if (!sensitivity_3) {
	if (staged_man.mortality_age_screen < staged_man.mortality_age_clinical)
	    staged_man.mortality_age_screen = staged_man.mortality_age_clinical;
    }

    return staged_man.mortality_age_screen - staged_man.mortality_age_clinical;
  };

  void apply_screen(stage_man& staged_man,
		    std::vector<int>& screen_ages,
		    Rcpp::List& cancer_mort_probs,
		    int staged_man_index,
		    int sensitivity_1)
  {
    int found_cancer_1 = 0;
    int found_cancer_2 = 0;
    int stage = -1;
    int stage_1 = -1;
    int stage_2 = -1;
    double screen_positive_age_1 = -1;
    double screen_positive_age_2 = -1;
    int first_screen = 0;
    for (int i: screen_ages) {
      stage = -1;
      int need_cancer_screen = 0;
      if (staged_man.diagnosis_age < i) {
	found_cancer_1 = 1;
	found_cancer_2 = 1;
      }
      if ((found_cancer_1 == 0) &&
	  staged_man.age_in_2021 <= i + 15 &&
	  staged_man.age_in_2021 > 44) {
        if (staged_man.screen_eligible) {
      canc_screens++;
	    screen_test_1(staged_man,
			i,
			need_cancer_screen);
	    screen_test_2(staged_man,
			i,
			need_cancer_screen);
	    if (!staged_man.positive_on_e6 && !staged_man.positive_on_hpv) {
	    number_of_infection_screens[0] += staged_man.weight;
	    } else if (!staged_man.positive_on_e6 && staged_man.positive_on_hpv) {
	    number_of_infection_screens[1] += staged_man.weight;
	    } else if (staged_man.positive_on_e6  && staged_man.positive_on_hpv) {
	    number_of_infection_screens[3] += staged_man.weight;
	    } else if (staged_man.positive_on_e6 && !staged_man.positive_on_hpv) {
	    number_of_infection_screens[2] += staged_man.weight;
	    }

	}

	if (!sensitivity_1) {
	    if (staged_man.positive_on_e6 || !staged_man.positive_on_hpv)
		staged_man.screen_eligible = 0;
	} else {
	  if (staged_man.positive_on_e6 || !staged_man.positive_on_hpv)
	    staged_man.screen_eligible = 0;
	  
	  if (!staged_man.positive_on_e6) { 
	    if (!staged_man.positive_on_hpv) {
		staged_man.negative_on_hpv++;
	    } else {
		staged_man.negative_on_hpv = 0;
	    }

	    if (first_screen > 0 && staged_man.negative_on_hpv < 2)
		staged_man.screen_eligible = 1;
	    }
	}

	
	first_screen++;
	// if (!staged_man.positive_on_navdx) {
	//     if (!staged_man.positive_on_e6 && !staged_man.positive_on_hpv) {
	//     number_of_infection_screens[4] += staged_man.weight;
	//     } else if (!staged_man.positive_on_e6 && staged_man.positive_on_hpv) {
	//     number_of_infection_screens[5] += staged_man.weight;
	//     } else if (staged_man.positive_on_e6  && staged_man.positive_on_hpv) {
	//     number_of_infection_screens[7] += staged_man.weight;
	//     } else if (staged_man.positive_on_e6 && !staged_man.positive_on_hpv) {
	//     number_of_infection_screens[6] += staged_man.weight;
	//     }

	// }

	int ultrasound_1 = 0;
	int ultrasound_2 = 0;
        if (staged_man.positive_on_e6 || staged_man.positive_on_hpv) {
	  if (!staged_man.positive_on_navdx){
	    screen_test_3(staged_man,
			  i,
			  need_cancer_screen);
	    }
	  number_of_cancer_screens[0] += staged_man.weight;
	    ultrasound_1 = 1;
          }

	  if (staged_man.positive_on_navdx) {
	    number_of_cancer_screens[1] += staged_man.weight;
	    ultrasound_2 = 1;
	  }

	  if (ultrasound_1 || ultrasound_2) {
	    for (int j = 0; j < 4; j++) {
		if (staged_man.stage_ages[j] <= i && staged_man.stage_ages[j] != 0) {
		    stage = j;
		}
	  }

	    double s2_scale = 1;
	    if (sensitivity_2) {
	      s2_scale = 0.5;
	    }
	    double r = R::runif(0.0, 1.0);
		if (stage > -1 && ((staged_man.N[stage] == 0 && staged_man.tumor_size[stage] == 0 && r < 0.0) ||
				   (staged_man.N[stage] == 0 && staged_man.tumor_size[stage] == 1 && r < s2_scale * 0.65) ||
				   (staged_man.N[stage] == 0 && staged_man.tumor_size[stage] == 2 && r < s2_scale * 0.75) ||
				   (staged_man.N[stage] >= 1 && r < s2_scale * 0.9))) {
		if (ultrasound_1 && found_cancer_1 == 0) {
		    number_of_cancer_screens[2] += staged_man.weight;
		    screen_positive_age_1 = i;
		    found_cancer_1 = 1;
		    stage_1 = stage;
		}

		if (ultrasound_2 && found_cancer_2 == 0) {
		    number_of_cancer_screens[3] += staged_man.weight;
		    screen_positive_age_2 = i;
		    found_cancer_2 = 1;
		    stage_2 = stage;
		}

	    }
        }
      }
    }

    int screen_1_age_diff;
    if (screen_positive_age_1 >= 0 && stage_1 > -1) {
      screen_age_stage_ajcc7(screen_positive_age_1, stage_1) += staged_man.weight;
      screen_age_stage_ajcc8(screen_positive_age_1, staged_man.ajcc8_stage[stage_1]) += staged_man.weight;
      if ((stage_1 < 2) ||
	  (stage_1 == 2 && staged_man.node_size[2] == 0 && R::runif(0.0, 1.0) < 0.67)) {
	tx_by_age(screen_positive_age_1, 0) += staged_man.weight;
	life_results_1(staged_man_index, 8) = 1;
      } else {
	tx_by_age(screen_positive_age_1, 1) += staged_man.weight;
	life_results_1(staged_man_index, 8) = 0;
      }

      for (int i = 0; i < screen_positive_age_1 + 1; i++)
	cancer_rate_denom[i] += staged_man.weight;
      cancer_rate_num[screen_positive_age_1] += staged_man.weight;

      screen_1_age_diff = get_mort_age_diff(screen_positive_age_1, staged_man,
					    Rcpp::as<Rcpp::List>(cancer_mort_probs[staged_man.smk])[stage_1],
					    Rcpp::as<Rcpp::List>(cancer_mort_probs[staged_man.smk])[staged_man.diagnostic_stage]);
    } else {
      screen_age_stage_ajcc7(staged_man.diagnosis_age, staged_man.diagnostic_stage) += staged_man.weight;
      screen_age_stage_ajcc8(staged_man.diagnosis_age, staged_man.ajcc8_stage[staged_man.diagnostic_stage]) += staged_man.weight;
      if ((staged_man.diagnostic_stage < 2) ||
	  (staged_man.diagnostic_stage == 2 && staged_man.node_size[2] == 0 && R::runif(0.0, 1.0) < 0.67)) {
	tx_by_age(staged_man.diagnosis_age, 0) += staged_man.weight;
	life_results_1(staged_man_index, 8) = 1;
      } else {
	tx_by_age(staged_man.diagnosis_age, 1) += staged_man.weight;
	life_results_1(staged_man_index, 8) = 0;
      }

      for (int i = 0; i < staged_man.diagnosis_age + 1; i++)
        cancer_rate_denom[i] += staged_man.weight;
      cancer_rate_num[staged_man.diagnosis_age] += staged_man.weight;

      stage_1 = staged_man.diagnostic_stage;
      screen_positive_age_1 = staged_man.diagnosis_age;
      screen_1_age_diff = get_mort_age_diff(screen_positive_age_1, staged_man,
					    Rcpp::as<Rcpp::List>(cancer_mort_probs[staged_man.smk])[stage_1],
					    Rcpp::as<Rcpp::List>(cancer_mort_probs[staged_man.smk])[staged_man.diagnostic_stage]);


    }
    if (screen_positive_age_2 >= 0 && stage_2 > -1) {
      screen_age_stage_ajcc7(screen_positive_age_2, stage_2 + 4) += staged_man.weight;
      screen_age_stage_ajcc8(screen_positive_age_2, staged_man.ajcc8_stage[stage_2] + 4) += staged_man.weight;
      if ((stage_2 < 2) ||
	  (stage_2 == 2 && staged_man.node_size[2] == 0 && R::runif(0.0, 1.0) < 0.67)) {
	tx_by_age(screen_positive_age_2, 2) += staged_man.weight;
      } else {
	tx_by_age(screen_positive_age_2, 3) += staged_man.weight;
      }


    } else {
      screen_age_stage_ajcc7(staged_man.diagnosis_age, staged_man.diagnostic_stage + 4) += staged_man.weight;
      screen_age_stage_ajcc8(staged_man.diagnosis_age, staged_man.ajcc8_stage[staged_man.diagnostic_stage] + 4) += staged_man.weight;
      if ((staged_man.diagnostic_stage < 2) ||
	  (staged_man.diagnostic_stage == 2 && staged_man.node_size[2] == 0 && R::runif(0.0, 1.0) < 0.67)) {
	tx_by_age(staged_man.diagnosis_age, 2) += staged_man.weight;
      } else {
	tx_by_age(staged_man.diagnosis_age, 3) += staged_man.weight;
      }

     }

    // difference in life expectancy (screen vs clinical)
    life_results_1(staged_man_index, 0) = screen_1_age_diff;
    // stage at detection (screen [if screen detected] otherwise clinical)
    life_results_1(staged_man_index, 1) = stage_1;
    // weight
    life_results_1(staged_man_index, 2) = staged_man.weight;
    // diference in age at diagnoses (clinical - screen)
    life_results_1(staged_man_index, 3) = staged_man.diagnosis_age - screen_positive_age_1;
    // age at clinical diagnosis
    life_results_1(staged_man_index, 4) = staged_man.diagnosis_age;
    // age at death in the absence of screening
    life_results_1(staged_man_index, 5) = staged_man.mortality_age_clinical;
    // age at death with screening
    life_results_1(staged_man_index, 6) = staged_man.mortality_age_screen;
    // clinical stage
    life_results_1(staged_man_index, 7) = staged_man.diagnostic_stage;
    // age in 2021
    life_results_1(staged_man_index, 9) = staged_man.age_in_2021;
    // screen death (from cancer)
    life_results_1(staged_man_index, 10) = staged_man.screen_death;
    // clinical death (from cancer)
    life_results_1(staged_man_index, 11) = staged_man.clinical_death;
    // age at all cause mortality
    life_results_1(staged_man_index, 12) = staged_man.mortality_age_no_cancer;
  };

  void apply_screen(infected_man& infection_man,
		    std::vector<int>& screen_ages,
		    int sensitivity_1)
  {
    int first_screen = 0;
    for (int i: screen_ages) {
      int need_cancer_screen = 0;
      if (infection_man.mortality > i && 
	  infection_man.age_in_2021 <= i + 15 &&
	  infection_man.age_in_2021 > 44) {

	if (infection_man.screen_eligible) {
	    inf_screens++;
	    screen_test_1(infection_man,
			i,
			need_cancer_screen);
	    screen_test_2(infection_man, i,
			need_cancer_screen);
	    first_screen++;
	    if (!infection_man.positive_on_e6 && !infection_man.positive_on_hpv) {
	    number_of_infection_screens[0] += infection_man.weight;
	    } else if (!infection_man.positive_on_e6 && infection_man.positive_on_hpv) {
	    number_of_infection_screens[1] += infection_man.weight;
	    } else if (infection_man.positive_on_e6  && infection_man.positive_on_hpv) {
	    number_of_infection_screens[3] += infection_man.weight;
	    } else if (infection_man.positive_on_e6 && !infection_man.positive_on_hpv) {
	    number_of_infection_screens[2] += infection_man.weight;
	    }

	}
	if (!sensitivity_1) {
	    if (infection_man.positive_on_e6 || !infection_man.positive_on_hpv)
		infection_man.screen_eligible = 0;
	} else {
	  if (infection_man.positive_on_e6 || !infection_man.positive_on_hpv)
	    infection_man.screen_eligible = 0;
	  
	  if (!infection_man.positive_on_e6) { 
	    if (!infection_man.positive_on_hpv) {
		infection_man.negative_on_hpv++;
	    } else {
		infection_man.negative_on_hpv = 0;
	    }

	    if (first_screen > 0 && infection_man.negative_on_hpv < 2)
		infection_man.screen_eligible = 1;
	    }
	}


	if (infection_man.positive_on_hpv || infection_man.positive_on_e6) {
	  number_of_cancer_screens[0] += infection_man.weight;
	}
      }
    }
  };

  void get_screen_ages(int age_in_2021,
		       std::vector<int>& screen_params,
		       std::vector<int>& screen_ages)
  {
    int current_age;
    if (age_in_2021 <= screen_params[0]) {
      current_age = screen_params[0] - 15;
    } else {
      current_age = age_in_2021 - 15;
    }


    int j = 1;
    if (current_age < 65) {
	if (screen_params[1] == 0) {
	    screen_ages.push_back(current_age);
	    if(screen_params[2] > 0) {
	      for (int i = 1; i < screen_params[2] + 1; i++) {
		screen_ages.push_back(current_age + i);
	      }
	    }
	} else {
	    while(current_age < 65) {
		screen_ages.push_back(current_age);
		current_age += screen_params[1];
		j++;
	    }
	}
    }

  };
  
  void sim_screen(std::vector<stage_man>& stage_data,
		  std::vector<infected_man>& infection_man_vec,
		  std::vector<int> screen_params,
		  Rcpp::List& cancer_mort_probs,
		  int sensitivity_1)
  {
    int staged_man_index = 0;
    for (auto iter = stage_data.begin(); iter != stage_data.end(); ++iter) {
      std::vector<int> screen_ages;
      if (screen_params.size() > 0) {
	get_screen_ages((*iter).age_in_2021,
			screen_params,
			screen_ages);
      }

      apply_screen((*iter), screen_ages, cancer_mort_probs, staged_man_index, sensitivity_1);
	screen_ages.clear();
	staged_man_index++;
    }

    for (auto iter = infection_man_vec.begin(); iter != infection_man_vec.end(); ++iter) {
      std::vector<int> screen_ages;
      if (screen_params.size() > 0) {
	get_screen_ages((*iter).age_in_2021,
			screen_params,
			screen_ages);
      }
      apply_screen((*iter), screen_ages, sensitivity_1);
      screen_ages.clear();
    }
  };

  void update_rates(const Rcpp::NumericMatrix& no_infection_mortality)
  {
    for (int i = 0; i < 70; i++)
      cancer_rate_denom[i] += no_infection_mortality(i, i);

    for (int i = 0; i < 70; i++)
      cancer_rate_num[i] /= cancer_rate_denom[i];
  }

  void write(int progression_speed,
	     int seed,
	     int scenario,
	     int iter_num)
  {
    std::string screen_file = "Results/cancer/screening_data/stage8/stage_age_" + std::to_string(progression_speed) + "_" +
      std::to_string(seed) + "_" + std::to_string(scenario) + "_" + std::to_string(iter_num) + ".csv";
    write_matrix(screen_file, screen_age_stage_ajcc8);
    
    std::string screen7_file = "Results/cancer/screening_data/stage7/stage_age_" + std::to_string(progression_speed) + "_" +
      std::to_string(seed) + "_" + std::to_string(scenario) + "_" + std::to_string(iter_num) + ".csv";
    write_matrix(screen7_file, screen_age_stage_ajcc7);

    std::string tx_file = "Results/cancer/screening_data/tx/tx_age_" + std::to_string(progression_speed) + "_" +
      std::to_string(seed) + "_" + std::to_string(scenario) + "_" + std::to_string(iter_num) + ".csv";
    write_matrix(tx_file, tx_by_age);

    Rcpp::NumericMatrix mat(2, 8);
    for (int i = 0; i < 8; i++)
	mat(0, i) = number_of_infection_screens[i];
    for (int i = 0; i < 4; i++)
      mat(1, i) = number_of_cancer_screens[i];

    std::string screen_num_file = "Results/cancer/screening_data/Nscreens/number_of_screens_" + std::to_string(progression_speed) + "_" +
      std::to_string(seed) + "_" + std::to_string(scenario) + "_" + std::to_string(iter_num) + ".csv";
    write_matrix(screen_num_file, mat);
    Rcpp::NumericMatrix cancer_rates(2, 70);
    for (int i = 0; i < 70; i++)
      cancer_rates(0, i) = cancer_rate_num[i];

//    std::string rate_file = "Results/cancer/screening_data/rates/rates_by_age_" + std::to_string(progression_speed) + "_" +
//      std::to_string(seed) + "_" + std::to_string(scenario) + "_" + std::to_string(iter_num) + ".csv";
//    write_matrix(rate_file, cancer_rates);
    
    std::string life_file = "Results/cancer/screening_data/life_expectancy/life_expectancy_" + std::to_string(progression_speed) + "_" +
      std::to_string(seed) + "_" + std::to_string(scenario) + "_" + std::to_string(iter_num) + ".csv";
    write_matrix(life_file, life_results_1);

  };

private:
  void write_matrix(std::string file_name,
                      const NumericMatrix& matrix)
    {
      std::ofstream outfile;
      outfile.open(file_name, std::ofstream::trunc);
         for (int i = 0; i < matrix.nrow(); i++) {
             for (int j = 0; j < matrix.ncol(); j++) {
                 outfile << matrix(i, j);
                 if (j == matrix.ncol() - 1) {
                     outfile << '\n';
                 } else {
                     outfile << ',';
                 }
             }
         }
         outfile.close();
    };
 
  
};

// [[Rcpp::export]]
int cpp_sim(Rcpp::NumericMatrix stage_probs,
	    Rcpp::NumericMatrix progression_params,
	    Rcpp::NumericMatrix cancer_data,
	    Rcpp::NumericMatrix infection_data,
	    Rcpp::NumericMatrix weights,
	    Rcpp::NumericMatrix mortality_data,
	    Rcpp::List TNM,
	    int progression_speed,
	    int seed,
	    Rcpp::List screen_scenarios,
	    int internal_iters,
	    Rcpp::NumericMatrix no_infection_mortality,
	    Rcpp::List cancer_mort_probs,
	    int sensitivity_1,
	    int sensitivity_2,
	    int sensitivity_3)
{
  int n = cancer_data.nrow();
  int m = infection_data.nrow();
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  
  std::vector<infected_man> infection_man_vec;
  set_seed_r(seed);
  for (int i = 0; i < m; i++) {
    infection_man_vec.push_back(infected_man(infection_data(i, _ ),
					   weights));
  }
  // no-screening
  std::vector<screen_data> screen_data_vec;
  for (int i = 0; i < internal_iters; i++) {
    screen_data_vec.push_back(screen_data(cancer_data.nrow(),
					  sensitivity_2,
					  sensitivity_3));
  }
  std::vector<int> no_screen{};
  int j = 0;
  for (auto iter = screen_data_vec.begin(); iter != screen_data_vec.end(); ++iter) {
  std::vector<stage_man> stage_data;
  set_seed_r(seed + j + 1);
  for (int i = 0; i < n; i++) {
  	stage_data.push_back(stage_man(cancer_data(i, _ ),
  					weights,
  					mortality_data,
  					stage_probs,
  					progression_params,
  				        progression_speed,
  				        TNM));
  }

  	(*iter).sim_screen(stage_data,
  			infection_man_vec,
			   no_screen,
			   cancer_mort_probs,
			   sensitivity_1);

	(*iter).update_rates(no_infection_mortality);
  	(*iter).write(progression_speed,
  		    seed,
		    0,
  		    j);
  	j++;
  	stage_data.clear();
    }
  screen_data_vec.clear();

  for (int k = 0; k < screen_scenarios.size(); k++) {
    std::vector<screen_data> screen_data_vec;
    for (int i = 0; i < internal_iters; i++) {
      screen_data_vec.push_back(screen_data(cancer_data.nrow(),
					    sensitivity_2,
					    sensitivity_3));
    }
    int j = 0;
    for (auto iter = screen_data_vec.begin(); iter != screen_data_vec.end(); ++iter) {
    std::vector<stage_man> stage_data;
    infection_man_vec.clear();
    set_seed_r(seed);
    for (int i = 0; i < m; i++) {
	infection_man_vec.push_back(infected_man(infection_data(i, _ ),
					    weights));
    }

    set_seed_r(seed + j + 1);
    for (int i = 0; i < n; i++) {
	stage_data.push_back(stage_man(cancer_data(i, _ ),
					weights,
					mortality_data,
					stage_probs,
					progression_params,
				        progression_speed,
				        TNM));
    }

	(*iter).sim_screen(stage_data,
			infection_man_vec,
			   screen_scenarios[k],
			   cancer_mort_probs,
			   sensitivity_1);

	(*iter).update_rates(no_infection_mortality);

	(*iter).write(progression_speed,
		    seed,
		    k + 1,
		    j);
	j++;

	stage_data.clear();
	}
    screen_data_vec.clear();
    
  }
  
  return 0;
}
