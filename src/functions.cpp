#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <string>

// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace std;

// Structure to hold data for simulated individual in prevalence sims
// We store all partner ages (in 8 categories) which do not change across cohorts
// We store a single set if incidence values which is updated with cohort
struct sim_man {
    NumericMatrix partner_ages = NumericMatrix(70, 8);
    NumericVector exposed = NumericVector(70);
    int first_age_of_partner;
    NumericVector new_partner_by_age = NumericVector(70);
    int has_cancer = 0;
    int cancer_age = -99;
    int cancer_infection_age;
    int ever_infected;
    int subcohort;
    int age_in_2021;
    int last_followup;
  double r_test;
  NumericVector e6_infections = NumericVector(70);
};

// Structure to hold data for simulated individual in cancer sims
// Now storing infection history
struct sim_man_cancer {
    sim_man sim_man;
    NumericVector infection_history = NumericVector(70);
    NumericVector infection_history_screen = NumericVector(70);
    NumericVector incident_infections = NumericVector(70);
  int ind;
};

// Structure holding all simulated data for single birth cohort
struct sim_birth_cohort {
    int internal_sim_N;
    NumericMatrix vaccine_results = NumericMatrix(internal_sim_N, 70);
    NumericMatrix cancer_results_no_smk = NumericMatrix(internal_sim_N, 70);
    NumericMatrix cancer_results_smk = NumericMatrix(internal_sim_N, 70);
    NumericMatrix cancer_after_vaccine_results = NumericMatrix(internal_sim_N, 70);
    NumericMatrix cancer_no_vaccine = NumericMatrix(internal_sim_N, 70);
    NumericMatrix hpv_prevalence = NumericMatrix(internal_sim_N, 70);
    NumericMatrix cancer_incidence = NumericMatrix(internal_sim_N, 70);
    NumericMatrix cancer_incidence_vaccine = NumericMatrix(internal_sim_N, 70);
    NumericMatrix cancer_incidence_no_vaccine = NumericMatrix(internal_sim_N, 70);
    NumericMatrix at_risk = NumericMatrix(internal_sim_N, 70);
    NumericMatrix vaccinated_at_risk = NumericMatrix(internal_sim_N, 70);
    NumericMatrix no_vaccinated_at_risk = NumericMatrix(internal_sim_N, 70);
    NumericMatrix hpv_incidence = NumericMatrix(internal_sim_N, 70);
    NumericMatrix hpv_incidence_vaccine = NumericMatrix(internal_sim_N, 70);
    NumericMatrix hpv_incidence_no_vaccine = NumericMatrix(internal_sim_N, 70);
    NumericMatrix hpv_prevalence_vaccine = NumericMatrix(internal_sim_N, 70);
    NumericMatrix hpv_prevalence_no_vaccine = NumericMatrix(internal_sim_N, 70);
    NumericMatrix incidence_at_risk = NumericMatrix(internal_sim_N, 70);
    NumericMatrix incidence_at_risk_vaccine = NumericMatrix(internal_sim_N, 70);
    NumericMatrix incidence_at_risk_no_vaccine = NumericMatrix(internal_sim_N, 70);


    sim_birth_cohort(int N) :
        internal_sim_N(N) {}
};

/* Calculate summary statistics of simulations */
/* summaries are basic means across internal simulations */
void get_summaries(const vector<sim_birth_cohort>& sim_data,
                   NumericMatrix& vaccine_means,
                   NumericMatrix& cancer_means,
                   NumericMatrix& cancer_vaccine_means,
                   NumericMatrix& cancer_no_vaccine_means,
                   NumericMatrix& cancer_incidence_means, 
                   NumericMatrix& cancer_incidence_vaccine_means, 
                   NumericMatrix& cancer_incidence_no_vaccine_means, 
                   NumericMatrix& at_risk_means, 
                   NumericMatrix& vaccinated_at_risk_means, 
                   NumericMatrix& no_vaccinated_at_risk_means, 
                   NumericMatrix& hpv_incidence_means, 
                   NumericMatrix& hpv_incidence_vaccine_means, 
                   NumericMatrix& hpv_incidence_no_vaccine_means, 
                   NumericMatrix& hpv_prevalence_means, 
                   NumericMatrix& hpv_prevalence_vaccine_means, 
                   NumericMatrix& hpv_prevalence_no_vaccine_means, 
                   NumericMatrix& incidence_at_risk_means, 
                   NumericMatrix& incidence_at_risk_vaccine_means, 
                   NumericMatrix& incidence_at_risk_no_vaccine_means, 
                   int internal_sim_N)
{
    for (int i = 0; i < 70; i++) {
	NumericMatrix cancer_results(internal_sim_N, 70);
	for (int k = 0; k < 70; k++) {
	    for (int l = 0; l < 70; l++) {
	    cancer_results(k, l) = sim_data[i].cancer_results_smk(k, l) + sim_data[i].cancer_results_no_smk(k, l);
	    }
	}
        for (int j = 0; j < 70; j++) {
            vaccine_means(i, j) = mean(sim_data[i].vaccine_results(_, j));
            cancer_means(i, j) = mean(cancer_results(_, j));
            cancer_vaccine_means(i, j) = mean(sim_data[i].cancer_after_vaccine_results(_, j));
            cancer_no_vaccine_means(i, j) = mean(sim_data[i].cancer_no_vaccine(_, j));
            cancer_incidence_means(i, j) = mean(sim_data[i].cancer_incidence(_, j));
            cancer_incidence_vaccine_means(i, j) = mean(sim_data[i].cancer_incidence_vaccine(_, j));
            cancer_incidence_no_vaccine_means(i, j) = mean(sim_data[i].cancer_incidence_no_vaccine(_, j));
            hpv_prevalence_means(i, j) = mean(sim_data[i].hpv_prevalence(_, j));
            hpv_incidence_means(i, j) = mean(sim_data[i].hpv_incidence(_, j));
            at_risk_means(i, j) = mean(sim_data[i].at_risk(_, j));
            vaccinated_at_risk_means(i, j) = mean(sim_data[i].vaccinated_at_risk(_, j));
            no_vaccinated_at_risk_means(i, j) = mean(sim_data[i].no_vaccinated_at_risk(_, j));
            hpv_incidence_vaccine_means(i, j) = mean(sim_data[i].hpv_incidence_vaccine(_, j));
            hpv_incidence_no_vaccine_means(i, j) = mean(sim_data[i].hpv_incidence_no_vaccine(_, j));
            hpv_prevalence_no_vaccine_means(i, j) = mean(sim_data[i].hpv_prevalence_no_vaccine(_, j));
            hpv_prevalence_vaccine_means(i, j) = mean(sim_data[i].hpv_prevalence_vaccine(_, j));
        }
    }
}

/* Calculate prevalence */
NumericVector get_prevalence(const NumericMatrix& states,
                             int n)
{  
    NumericVector prev(70);
    
    for (int i = 0; i < 70; i++)
        prev[i] = sum(states(i, _)) / n;
    
    return prev;
}

/* Calculate cumulative incidence estimates */
NumericVector get_cum_inc(const List& cum_inc_list,
                       const NumericMatrix& weights,
                       const double n)
{
    NumericVector raw_cum_inc(70);

    for (int i = 0; i < 32; i++) {
        NumericVector cum_inc_tmp = cum_inc_list[i];
        for (int j = 0; j < 70; j++) {
            raw_cum_inc[j] += weights(j, i) * cum_inc_tmp[j] / n;
        }

    }

    return raw_cum_inc;
}

/* Calculate logistic probability for clearance calculation */
double pweibull(double x, double shape, double scale)
{
    if (x < 0)
        return 0;

    return R::pweibull(x, shape, scale, 1, 0);
}

// Calculate clearance probability at point x
// Values shifted by 0.5 years for mid year acquisition
double gen_p_clear(double x, double shape, double scale) 
{

  return (pweibull((x+0.5), shape, scale) - pweibull( (x - 0.5), shape, scale)) /
    (1 - pweibull((x - 0.5), shape, scale));
}

/* Generate vector of probabilities of clearance by duration of infection */
NumericVector get_p_clear(double shape,
                          double scale)
{
    NumericVector p_clear(70);

    for (int i = 0; i < 70; i++) {
        p_clear[i] = gen_p_clear(i, shape, scale);
    }

    return p_clear;
}

/* Generate multinomial distribution with uniform probabilities.
 * Used to assign specific partner age. */
inline int gen_multi_int(int min, int max) 
{
     int range = max - min + 1;
     double r = R::runif(0.0, 1.0);
     int val = 0;
     for (int j = 1; j < range; j++) {
          if (r > (double) j / range) {
               val++;
          } else {
               break;
          }
     }
     val += min;

     return val;
}

/* Get specific partner age from grouped variable */
inline int get_part_age(int male_age, int part_group)
{
     int female_age = 0;

     if (part_group == 0) {
          female_age = male_age - gen_multi_int(15, 20);
     } else if (part_group == 1) {
          female_age = male_age - gen_multi_int(10, 14);
     } else if (part_group == 2) {
          female_age = male_age - gen_multi_int(5, 9);
     } else if (part_group == 3) {
          female_age = male_age - gen_multi_int(0, 4);
     } else if (part_group == 4) {
          female_age = male_age + gen_multi_int(1, 4);
     } else if (part_group == 5) {
          female_age = male_age + gen_multi_int(5, 9);
     } else if (part_group == 6) {
          female_age = male_age + gen_multi_int(10, 14);
     } else if (part_group == 7) {
          female_age = male_age + gen_multi_int(15, 20);
     }

     // Note: For random assignment with groups under 15, currently weighted to
     // 15.
     if (female_age < 0) {
          female_age = 0;

    }

     return female_age;
}

// Generate incidence values for individual at each age
void get_exposed(sim_man& sim_man,
                           const NumericVector& exposure_age_dist,
                           const NumericMatrix& herd_immunity,
		 int start_year)
{
    int exposure;
    double herd_immunity_value;
    int cur_year = start_year;
   
    for (int j = 0; j < 70; j++) {
        exposure = 0;
        for (int i = 0; i < 8; i++) {
            int l = sim_man.partner_ages(j, i);

            if (l == 0)
                continue;

            while (l > 0) {
                int part_age = get_part_age(j, i);

                if (part_age <= 45) {
                    herd_immunity_value = herd_immunity(cur_year, part_age);
                } else {
                    herd_immunity_value = herd_immunity(cur_year - (part_age - 45), 45);
                }

                double r = R::runif(0.0, 1.0);
                exposure += r < (1-herd_immunity_value) * 
                    exposure_age_dist[part_age],
                    l--;
            }
        }
        sim_man.exposed[j] = exposure;
        cur_year++;
    }
}

/* Assign age to each partner */
inline NumericVector get_part_age_vec(int male_age,
                                    const int& num_parts, 
                                    const NumericVector& partner_age_dist)
{
    double r;
    int val;
    NumericVector part_age_vec(8);

    for (int i = 0; i < num_parts; i++) {
        r = R::runif(0.0, 1.0);
        val = 0; 
        // First assign an age group relative to male
        for (int j = 0; j < 7; j++) {
            if (r > partner_age_dist[j]) {
                val++;
            }
        }

        part_age_vec[val]++;
    }

    return part_age_vec;
}

/* Draw number of partners */
inline int get_num_partners(const double& num_partner_dist_mn,
                     const NumericVector& num_partner_dist_pois,
                     const int age,
                     sim_man& sim_man)
{
    double r = R::runif(0.0, 1.0);
//    int num_parts = 0;
    NumericVector num_parts;

    if (r < num_partner_dist_mn) {
         sim_man.new_partner_by_age[age] = 1;
	// First and Third column: size; Second and Fourth column: mu
	num_parts = Rcpp::rnbinom_mu(1,num_partner_dist_pois[2], num_partner_dist_pois[3]);
    } else {
	num_parts = Rcpp::rnbinom_mu(1,num_partner_dist_pois[0], num_partner_dist_pois[1]);
    }

    return num_parts[0];
}

/* Generate number of partners by age */
void sim_partners(sim_man& sim_man,
                  const NumericVector& num_partner_dist_mn,
                  const NumericMatrix& num_partner_dist_pois,
                  const NumericMatrix& partner_age_dist)
{
    int num_partners;
    int has_a_partner = 0;
    for (int i = 0; i < 70; i++) {
        num_partners = get_num_partners(num_partner_dist_mn[i],
                                        num_partner_dist_pois(i, _),
                                        i,
                                        sim_man);
        if (has_a_partner == 0) {
             if (num_partners > 0) {
                  sim_man.first_age_of_partner = i;
                  has_a_partner++;
             }
        }
        sim_man.partner_ages(i, _) = get_part_age_vec(i, num_partners, 
                                                      partner_age_dist(i, _));
    }

}

// For natural history simulations, simulate full natural history for
// single individual
void sim_states(sim_man& sim_data,
                const int cohort,
                const NumericMatrix& weights,
                NumericMatrix& states,
                const NumericVector& p_clear,
                const NumericVector& exposure_age_dist,
                double trans_prob,
                double immune_scale,
                NumericVector& cum_inc,
                NumericVector& incidence)
{
    int infected = 0;
    bool ever_infected = 0;
    double trans_prob_tmp = trans_prob;

    for (int i = 0; i < 70; i++) {
        if (ever_infected)
             trans_prob_tmp = trans_prob * immune_scale;
            //sim_data.incidence[i] *= immune_scale;
        
        if (infected > 0) {
             if (R::runif(0.0, 1.0) < p_clear(infected)) {
                 infected = 0;
             } else {
                 infected += 1;
             }
        }

        if (infected == 0 && sim_data.exposed[i] > 0) {
             for (int j = 1; j < sim_data.exposed[i] + 1; j++) {
                 infected = R::runif(0.0, 1.0) < trans_prob_tmp;
                 if (infected) {
                      incidence[i] += weights(i, cohort);
                      ever_infected = 1;
                      if (R::runif(0.0, 1.0) < p_clear(0)) {
                           infected = 0;
                      } else {
                           break;
                      }
                 }
             }
        } 

        if (infected > 0) {
            states(i, infected - 1) += weights(i, cohort);
        }

        if (ever_infected) {
            cum_inc[i]++;
        }

	sim_data.e6_infections[i] = infected;
    }
}

/* For cancer simulations, simulate hpv history for single individual */
void sim_hpv(sim_man_cancer& sim_data,
                const NumericVector& p_clear,
                double trans_prob,
                double immune_scale)
{

    int infected = 0;
    int ever_infected = 0;
    double trans_prob_tmp = trans_prob;

    for (int i = 0; i < 70; i++) {
       if (ever_infected)
             trans_prob_tmp = trans_prob * immune_scale;
            //sim_data.incidence[i] *= immune_scale;
        
        if (infected > 0) {
             if (R::runif(0.0, 1.0) < p_clear(infected)) {
                 infected = 0;
             } else {
                 infected += 1;
             }
        }

        if (infected == 0 && sim_data.sim_man.exposed[i] > 0) {
             for (int j = 1; j < sim_data.sim_man.exposed[i] + 1; j++) {
                 infected = R::runif(0.0, 1.0) < trans_prob_tmp;
                 if (infected) {
                      sim_data.incident_infections(i)++;
                      ever_infected = 1;
                      if (R::runif(0.0, 1.0) < p_clear(0)) {
                           infected = 0;
                      } else {
                           break;
                      }
                 }
             }
        }

       sim_data.infection_history(i) = infected;
    }

}

/* Simulate natural histories for cancer simulations.*/
void sim_states(sim_man_cancer& sim_man,
                sim_birth_cohort& sim_cohort_data,
                const NumericVector& penetrance,
                const NumericMatrix& herd_immunity,
                const NumericMatrix& mortality_matrix,
                const NumericVector& vaccine_probs,
                const NumericMatrix& weights,
                int j, // Internal sim number
                int t, // Birth cohort indicator
                int cohort,
                NumericMatrix& subcohort_incidence,
		NumericMatrix& mort_denom,
		const int smk)
{
     int mortality = 0;
     int vaccine = 0;
     int vaccine_effective = 0;
     NumericVector infection_history(70);
     NumericVector incident_infections(70);
     for (int i = 0; i < 70; i++) {
          infection_history[i] = sim_man.infection_history[i];
          incident_infections[i] = sim_man.incident_infections[i];
     }
     sim_man.sim_man.has_cancer = 0;
     sim_man.sim_man.cancer_infection_age = 0;
     sim_man.sim_man.ever_infected = 0;
     int found_cancer = 0;
     int mortality_age = 69;

     for (int i = 0; i < 70; i++) {
          // Check if dead, otherwise include in inc/prev denominator
          if (R::runif(0.0, 1.0) < mortality_matrix(i, cohort) && (i + t > 68)) {
               mortality = 1;
	    mortality_age = i;
	    if (mortality_age > 69)
	      Rcout << i << " " << mortality_age << std::endl;
          }

          if (mortality) {
               break;
          } else {
               sim_cohort_data.at_risk(j, i) += weights(69 - t, cohort);
	  }

          // Check if vaccinated and vaccine efficacy
          if (!vaccine && R::runif(0.0, 1.0) < vaccine_probs[i]) {
               vaccine = 1;
               sim_cohort_data.vaccine_results(j, i) += weights(69 - t, cohort);
               if (R::runif(0.0, 1.0) < 0.95) {
                    vaccine_effective = 1;
                    for (int k = i; k < 70; k++)
                         incident_infections[k] = 0;
               }
          }

          if (vaccine_effective && (infection_history[i] < 2)) {
               for (int k = i; k < 70; k++)
                    infection_history[k] = 0;
          }

          if (vaccine) {
               sim_cohort_data.vaccinated_at_risk(j, i) += weights(69 - t, cohort);
          } else {
               sim_cohort_data.no_vaccinated_at_risk(j, i) += weights(69 - t, cohort);
          }

          // Increment incidence estimates if infections occured
          if (incident_infections[i] > 0) {
               sim_cohort_data.hpv_incidence(j, i) += incident_infections[i] *
                                                  weights(69 - t, cohort);
		    if (vaccine) {
                    sim_cohort_data.hpv_incidence_vaccine(j, i) += incident_infections[i] *
                                                  weights(69 - t, cohort);
               } else {
                    sim_cohort_data.hpv_incidence_no_vaccine(j, i) += incident_infections[i] *
                                                  weights(69 - t, cohort);
               }

		    if (i + t > 68) {
		      sim_man.sim_man.ever_infected = 1;
		    }
          }

          if (infection_history[i] < 2) {
               sim_cohort_data.incidence_at_risk(j, i) += weights(69 - t, cohort);

               if (vaccine) {
                    sim_cohort_data.incidence_at_risk_vaccine(j, i) += weights(69 - t, cohort);
               } else {
                    sim_cohort_data.incidence_at_risk_no_vaccine(j, i) += weights(69 - t, cohort);
               }

               if (infection_history[i] == 1 && (sim_man.sim_man.has_cancer == 0))
                    sim_man.sim_man.cancer_infection_age = i;
          }

          if (infection_history[i] > 0) {
               sim_cohort_data.hpv_prevalence(j, i) += weights(69 - t, cohort);
               
               if (vaccine) {
                    sim_cohort_data.hpv_prevalence_vaccine(j, i) += weights(69 - t, cohort);
               } else {
                    sim_cohort_data.hpv_prevalence_no_vaccine(j, i) += weights(69 - t, cohort);
               }

               if (R::runif(0.0, 1.0) > 1 - penetrance[infection_history[i] - 1]) {
		 if(smk) {
                    sim_cohort_data.cancer_results_smk(j, i) += weights(69 - t, cohort);
		 } else {
                    sim_cohort_data.cancer_results_no_smk(j, i) += weights(69 - t, cohort);
		 }
                    sim_man.sim_man.has_cancer = 1;
                    sim_man.sim_man.cancer_age = i;
                    mortality = 1;               
		    found_cancer = 1;
                    if (vaccine) {
                         sim_cohort_data.cancer_after_vaccine_results(j, i) += weights(69 - t, cohort);
                         sim_cohort_data.cancer_incidence_vaccine(j, i) += weights(69 - t, cohort);
                    } else {
                         sim_cohort_data.cancer_no_vaccine(j, i) += weights(69 - t, cohort);
                         sim_cohort_data.cancer_incidence_no_vaccine(j, i) += weights(69 - t, cohort);
                    }
               }
	       if (i + t == 69)
		 sim_man.sim_man.ever_infected = 1;
          }

	  if (infection_history[i] == 0 && incident_infections[i] > 0) {
	    sim_man.infection_history[i] = 0.5;
	  }

          // Subcohort incidence calculations
          if (vaccine == 0 && i < 31 && i > 11 && i + t > 68) {
               if (incident_infections[i] > 0) {
                    // No weights since looking at risk within subcohort
                    // Will need weights to make concentration curve
                    // (number of people in each subcohort)
                    subcohort_incidence(j, 0) += incident_infections[i];
               }

               if (infection_history[i] < 2) {
                    subcohort_incidence(j, 1)++;
               }
          }
          sim_man.sim_man.last_followup = i;
     }

     if (!found_cancer) {
	 for (int i = 0; i < mortality_age + 1; i++) {
	    mort_denom(t, i) += weights(69 - t, cohort);
	 }
     }
	 

     if (sim_man.sim_man.has_cancer == 0) {
          sim_man.sim_man.cancer_infection_age = -99;
          sim_man.sim_man.cancer_age = -99;
     }
}

/* For natural history simulations: Generate full simulation data
 * for single group.*/
void sim_cohort(vector<sim_man>& sim_data,
             const NumericVector& num_partner_dist_mn,
             const NumericMatrix& num_partner_dist_pois,
             const NumericMatrix& partner_age_dist,
             const int cohort,
             const NumericMatrix& weights,
             NumericMatrix& states,
             const NumericVector& p_clear,
             const NumericVector& exposure_age_dist,
             double trans_prob,
             double immune_scale,
             List& cum_inc_list,
             NumericVector& seed_vec,
		NumericVector& incidence,
		NumericMatrix& e6_results)
{
    NumericMatrix herd_immunity(70, 80);

    NumericVector cum_inc(70);
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];

    // Simulate partners/exposure seperately for random number consistency
    int t = 0;
    for (auto iter = sim_data.begin(); iter != sim_data.end(); ++iter) {
        set_seed_r(seed_vec[t] + cohort);
        sim_partners(*iter,
                     num_partner_dist_mn,
                     num_partner_dist_pois,
                     partner_age_dist);

        t++;
    }

    t = 0;
    for (auto iter = sim_data.begin(); iter != sim_data.end(); ++iter) {
        set_seed_r(seed_vec[t] + cohort);
        get_exposed(*iter,
                      exposure_age_dist,
                      herd_immunity,
		    0);
        sim_states(*iter,
                   cohort,
                   weights,
                   states,
                   p_clear,
                   exposure_age_dist,
                   trans_prob,
                   immune_scale,
                   cum_inc,
                   incidence);

	for (int i = 1; i < 70; i++) {
	  for (int j = 0; j < 70; j++) {
	    if ((*iter).e6_infections[j] == i) {
	      for (int k = j; k < 70; k++) {
		e6_results(k, i) += weights(k, cohort);
	      }
	      break;
	    }
	  }
	}

		
        t++;
       
    }


    cum_inc_list.push_back(cum_inc);
}

/* For cancer simulations: Generate full simulation data
 * for single group.*/
void sim_cohort(vector<sim_man_cancer>& sim_data,
             vector<sim_birth_cohort>& sim_cohort_data,
             const NumericVector& num_partner_dist_mn,
             const NumericMatrix& num_partner_dist_pois,
             const NumericMatrix& partner_age_dist,
             const int cohort,
             const NumericMatrix& weights,
             const NumericVector& p_clear,
             const NumericVector& exposure_age_dist,
             double trans_prob,
             double immune_scale,
             const NumericMatrix& herd_immunity,
             const NumericMatrix& vaccine_probs,
             const NumericVector& penetrance,
             const NumericMatrix& mortality,
             int internal_sim_N,
             NumericVector seed_vec,
             vector<sim_man_cancer>& screening_data,
		NumericMatrix& no_infections,
		NumericMatrix& no_infections_five_before,
		NumericMatrix& mort_denom,
		const int smk)
{

    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
  
    int t = 0;
    for (auto iter = sim_data.begin(); iter != sim_data.end(); ++iter) {
             set_seed_r(seed_vec[t] + cohort);
             sim_partners((*iter).sim_man,
                          num_partner_dist_mn,
                          num_partner_dist_pois,
                          partner_age_dist);

          t++; 
    }
    int s = 0; // birth cohort counter
   // First column: incident infections, second: at risk
   NumericMatrix subcohort_incidence = NumericMatrix(internal_sim_N, 2);
    for (auto iter_year = sim_cohort_data.begin(); iter_year != sim_cohort_data.end(); ++iter_year) {
        t = 0; // individual counter
        for (auto iter = sim_data.begin(); iter != sim_data.end(); ++iter) {
             set_seed_r(seed_vec[t] + cohort);
               get_exposed((*iter).sim_man,
                             exposure_age_dist,
                             herd_immunity,
			   s);
               sim_hpv(*iter,
                    p_clear,
                    trans_prob,
                    immune_scale);

            for (int j = 0; j < internal_sim_N; j++) {
                sim_states(*iter,
                        *iter_year,
                        penetrance,
                        herd_immunity,
                        mortality,
                        vaccine_probs(s, _),
                        weights,
                        j,
                        s,
                        cohort,
			   subcohort_incidence,
			   mort_denom,
			   smk);
            }

            if ((*iter).sim_man.ever_infected == 1) {

                 sim_man_cancer screen_man;
                 /* int index = (strat_data_index + cohort) * sim_data.size() + t; */
                 screen_man.sim_man.age_in_2021 = 69 - s + 15; // age in 2021
                 screen_man.sim_man.subcohort = cohort;
                 screen_man.sim_man.cancer_age = (*iter).sim_man.cancer_age;
                 screen_man.sim_man.cancer_infection_age = (*iter).sim_man.cancer_infection_age;
                 screen_man.sim_man.has_cancer = (*iter).sim_man.has_cancer;
                 screen_man.sim_man.subcohort = (*iter).sim_man.subcohort;
                 screen_man.sim_man.last_followup = (*iter).sim_man.last_followup;
		 screen_man.ind = t;
                screening_data.push_back(screen_man); 
                for (int j = 0; j < 70; j++){
                     screening_data[screening_data.size() - 1].infection_history_screen[j] = (*iter).infection_history[j];
                }

            } else {
	      int five_before = 0;
	    double r = R::rnorm(13, 5 / 1.96);
	    if (r < 1)
	      r = 1;

	      for (int k = 0; k < 69 - s + 1; k++) {
		// Generate year value from normal distribution. Min value of 1. In place of 13
		// Mean: 13, sd: 5 / 1.96
		if ((*iter).infection_history[k] > r) {
		  five_before = 1;
		  break;
		}
	      }
	      if (five_before) {
		no_infections_five_before(cohort, s)++;
	      } else {
                 no_infections(cohort, s)++;
	      }
            }

            for (int j = 0; j < 70; j++) {
               (*iter).incident_infections[j] = 0;
               (*iter).infection_history[j] = 0;
            }

            t++;
        }

	    Rcout << cohort << " : " << s << endl; 
        s++;
    }
}

/* Calculate rates for cancer history simulations. Number at risk is variable
 * for these calculations.*/
void calculate_rates(sim_birth_cohort sim_data,
                     int internal_sim_N)
{
	NumericMatrix cancer_results(internal_sim_N, 70);
	for (int k = 0; k < 70; k++) {
	    for (int l = 0; l < 70; l++) {
	    cancer_results(k, l) = sim_data.cancer_results_smk(k, l) + sim_data.cancer_results_no_smk(k, l);
	    }
	}
        for (int i = 0; i < internal_sim_N; i++) {
            for (int j = 0; j < 70; j++) {
                 if (sim_data.at_risk(i,j) > 0) {
                      if (sim_data.incidence_at_risk(i, j) > 0) {
                         sim_data.hpv_incidence(i, j) 
                              /= sim_data.incidence_at_risk(i, j);
                      }
                     sim_data.hpv_prevalence(i, j) 
                          /= sim_data.at_risk(i, j);
		    sim_data.cancer_incidence(i, j)
                          = cancer_results(i, j) / sim_data.at_risk(i, j);
                 } else {
                     sim_data.hpv_prevalence(i, j) 
                          = 0;
                     sim_data.cancer_incidence(i, j)
                          = 0;
                 }
                 if (sim_data.vaccinated_at_risk(i,j) > 0) {
                      if (sim_data.incidence_at_risk_vaccine(i, j) > 0) {
                         sim_data.hpv_incidence_vaccine(i, j) 
                              /= sim_data.incidence_at_risk_vaccine(i, j);
                      }
                     sim_data.hpv_prevalence_vaccine(i, j) 
                          /= sim_data.vaccinated_at_risk(i, j);
                     sim_data.cancer_incidence_vaccine(i, j)
                          = sim_data.cancer_incidence_vaccine(i, j) / sim_data.vaccinated_at_risk(i, j);
                 } else {
                     sim_data.hpv_prevalence_vaccine(i, j) 
                          = 0;
                     sim_data.cancer_incidence_vaccine(i, j)
                          = 0;
                 }                 
                 if (sim_data.no_vaccinated_at_risk(i,j) > 0) {
                      if (sim_data.incidence_at_risk_no_vaccine(i, j) > 0) {
                         sim_data.hpv_incidence_no_vaccine(i, j) 
                              /= sim_data.incidence_at_risk_no_vaccine(i, j);
                      }
                     sim_data.hpv_prevalence_no_vaccine(i, j) 
                          /= sim_data.no_vaccinated_at_risk(i, j);
                     sim_data.cancer_incidence_no_vaccine(i, j)
                          = sim_data.cancer_incidence_no_vaccine(i, j) / sim_data.no_vaccinated_at_risk(i, j);
                 } else {
                     sim_data.hpv_prevalence_no_vaccine(i, j) 
                          = 0;
                     sim_data.cancer_incidence_no_vaccine(i, j)
                          = 0;
                 }
            }
        }
}

/* Store natural history output as csv file when needed. */
void write_matrix(string file_name,
                  const NumericMatrix& matrix)
{
       ofstream outfile;
       outfile.open(file_name, ofstream::trunc);
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
}

// Primary function called for natural history simulations
// [[Rcpp::export]]
List sim_iter(int seed,
             NumericVector seed_vec,
             int n,
             List num_partner_dist_mn,
             List num_partner_dist_pois,
             NumericMatrix partner_age_dist,
             NumericMatrix weights,
             NumericMatrix exposure_age_dist,
             int write_states,
             double immune_scale,
             NumericVector clearance_scale,
	     int write_smoking)
{
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];

    // Create vector to store simulated data
    vector<sim_man> sim_data(n);

    // Create matrix to store weighted infection totals for non-smokers
    NumericMatrix states(70, 70);

    // Create matrix to store weighted infection totals for smokers
    NumericMatrix states_smk(70, 70);

    // Create vector to store incidence
    NumericVector incidence(70);

    // Initiate vector to store prevalence
    NumericVector prevalence(70);
   
    // Create vector to store cumulative incidence
    List cum_inc_list = List::create();
    NumericVector cum_inc(70);

    // set seed using first seed value
    set_seed_r(seed);
    
    double trans_prob = R::runif(0.0, 1.0);

    // Generate parameters for clearance distribution
    double shape = R::runif(0.1754,1.2280);
    double scale = R::runif(0.1719, 1.2032);

    // Generate clearance probabilities by length of infection
    // First value: 0-6 month clearance; Second: 6-18 months; Third: 18-30
    // months, etc.
    NumericVector p_clear = get_p_clear(shape, scale);

    // Additional clearance vector for smokers
    NumericVector p_clear_smk(p_clear.size());

	if(clearance_scale.size() != p_clear.size())
		Rcout << "Clearance scale is the wrong length!@!!!!" << endl;

    for (int i = 0; i < p_clear.size(); i++) {
	p_clear_smk[i] = clearance_scale[i] * p_clear[i];
    }

    NumericVector p_clear_tmp = p_clear;

    NumericMatrix e6_results(70, 70);

    // i: race
    // j: smoking status
    // k: quartile
     for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 4; k++) {
                if (j == 1) {
                      p_clear_tmp = p_clear_smk;

			sim_cohort(sim_data,
				   num_partner_dist_mn[8 * i + 4 * j + k],
				   num_partner_dist_pois[8 * i + 4 * j + k],
				   partner_age_dist,
				   8 * i + 4 * j + k, // Cohort value
				   weights,
				   states_smk,
				   p_clear_tmp,
				   exposure_age_dist(k, _),
				   trans_prob,
				   immune_scale,
				   cum_inc_list,
				   seed_vec,
				   incidence,
				   e6_results);
                 } else {
                      p_clear_tmp = p_clear;
			sim_cohort(sim_data,
				   num_partner_dist_mn[8 * i + 4 * j + k],
				   num_partner_dist_pois[8 * i + 4 * j + k],
				   partner_age_dist,
				   8 * i + 4 * j + k, // Cohort value
				   weights,
				   states,
				   p_clear_tmp,
				   exposure_age_dist(k, _),
				   trans_prob,
				   immune_scale,
				   cum_inc_list,
				   seed_vec,
				   incidence,
				   e6_results);
		}
            }
        }
    }

    for (int i = 0; i < 70; i++) {
      for (int j = 0; j < 70; j++) {
	e6_results(i, j) /= sim_data.size();
      }
    }

   cum_inc = get_cum_inc(cum_inc_list, weights, n);
    // Combine state matrices for prevalence
    NumericMatrix states_prev(70,70);
   for (int i = 0; i < 70; i++) {
	for (int j = 0; j < 70; j++) {
		states_prev(i, j) = states(i,j) + states_smk(i, j);
	}
   }

    prevalence = get_prevalence(states_prev, n);

    if(write_states) {
	if (write_smoking) {
        string file_name = "Results/states/states_no_smk" + to_string(seed) + ".csv";
        ofstream outfile;
        outfile.open(file_name, ofstream::trunc);
        for (int i = 0; i < 70; i++) {
            for (int j = 0; j < 70; j++) {
                outfile << states(i, j);
                if (j == 69) {
                    outfile << '\n';
                } else {
                    outfile << ',';
                }
            }
        }
        outfile.close();
        string file_name2 = "Results/states/states_smk" + to_string(seed) + ".csv";
        ofstream outfile2;
        outfile2.open(file_name2, ofstream::trunc);
        for (int i = 0; i < 70; i++) {
            for (int j = 0; j < 70; j++) {
                outfile2 << states_smk(i, j);
                if (j == 69) {
                    outfile2 << '\n';
                } else {
                    outfile2 << ',';
               }
            }
        }
        outfile2.close();
	} else {
		string file_name = "Results/states/states_" + to_string(seed) + ".csv";
		ofstream outfile;
		outfile.open(file_name, ofstream::trunc);
		for (int i = 0; i < 70; i++) {
		    for (int j = 0; j < 70; j++) {
			outfile << states_prev(i, j);
			if (j == 69) {
			    outfile << '\n';
			} else {
			    outfile << ',';
			}
		    }
		}
		outfile.close();
	}
    }

    string file_name = "Results/e6_check_" + to_string(seed) + ".csv";
    write_matrix(file_name, e6_results);
    

    return List::create(prevalence, shape, scale, trans_prob,
                        p_clear, cum_inc);
}

// Primary function called for cancer simulations
// [[Rcpp::export]]
int sim_iter_cancer_cpp(int vac_scenario,
             int seed,
             NumericVector seed_vec,
             int n,
             List num_partner_dist_mn,
             List num_partner_dist_pois,
             NumericMatrix partner_age_dist,
             NumericMatrix weights,
             NumericMatrix exposure_age_dist,
             double immune_scale,
             NumericVector penetrance_no_smk,
             NumericVector penetrance_smk,
             NumericMatrix herd_immunity,
             NumericMatrix vaccine_probs,
             NumericMatrix mortality,
             int internal_sim_N,
             NumericVector clearance_scale)
{
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    // Create vector to store simulated infection data
    vector<sim_man_cancer> sim_data(n);

    // Create vector to store cohort data
    vector<sim_birth_cohort> sim_cohort_data;
    for (int i = 0; i < 70; i++)
      sim_cohort_data.push_back(sim_birth_cohort(internal_sim_N));


    // set seed
    set_seed_r(seed);
 
    double trans_prob = R::runif(0.0, 1.0);

    // Generate parameters for clearance distribution
    double shape = R::runif(0.1754,1.2280);
    double scale = R::runif(0.1719, 1.2032);

    // Generate clearance probabilities by length of infection
    // First value: 0-6 month clearance; Second: 6-18 months; Third: 18-30
    // months, etc.
    NumericVector p_clear = get_p_clear(shape, scale);

    // Additional clearance vector for smokers
    NumericVector p_clear_smk(p_clear.size());

	if(clearance_scale.size() != p_clear.size())
		Rcout << "Clearance scale is the wrong length!@!!!!" << endl;

    for (int i = 0; i < p_clear.size(); i++) {
         p_clear_smk[i] = clearance_scale[i] * p_clear[i];
    }

    NumericVector p_clear_tmp = p_clear;    

    NumericVector penetrance_tmp = penetrance_no_smk;

    // Matrix to store data by subcohort
         // Matrix to store risk stratification data
    NumericMatrix no_infections(32,70);
    NumericMatrix no_infections_five_before(32, 70);

    vector<sim_man_cancer> screening_data;
    NumericMatrix mort_denom(70, 70);
    for (int i = 0; i < 4; i++) {
           for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 4; k++) {
		  screening_data.clear();

                  if (j == 1) {
                      p_clear_tmp = p_clear_smk;
                      penetrance_tmp = penetrance_smk;
                  } else {
                       p_clear_tmp = p_clear;
                       penetrance_tmp = penetrance_no_smk;
                  }
                   sim_cohort(sim_data,
                              sim_cohort_data,
                               num_partner_dist_mn[8 * i + 4 * j + k],
                               num_partner_dist_pois[8 * i + 4 * j + k],
                               partner_age_dist,
                               8 * i + 4 * j + k, //cohort value
                               weights,
                               p_clear_tmp,
                               exposure_age_dist(k, _),
                               trans_prob,
                               immune_scale,
                               herd_immunity,
                               vaccine_probs,
                               penetrance_tmp,
                               mortality,
                               internal_sim_N,
                               seed_vec,
                               screening_data,
			      no_infections,
			      no_infections_five_before,
			      mort_denom,
			      j);

	    int screen_rows = screening_data.size();

	NumericMatrix screen_matrix(screen_rows, 77);

	int t = 0;
	for (auto iter = screening_data.begin(); iter != screening_data.end(); ++iter) {
	    screen_matrix(t, 0) = 8 * i + 4 * j  + k;
	    screen_matrix(t, 1) = (*iter).sim_man.age_in_2021;
	    screen_matrix(t, 2) = (*iter).sim_man.has_cancer;
	    screen_matrix(t, 3) = (*iter).sim_man.cancer_infection_age;
	    screen_matrix(t, 4) = (*iter).sim_man.cancer_age;
	    screen_matrix(t, 5) = (*iter).sim_man.last_followup;

	    for (int i = 0; i < 70; i++) {
		screen_matrix(t, 6 + i) = (*iter).infection_history_screen[i];
	    }
	    screen_matrix(t, 76) = (*iter).ind;
	    t++;
	}
        string screen_file = "Results/cancer/screening_data" + to_string(vac_scenario) + "_" + to_string(seed) + "/" + to_string(8 * i + 4 * j + k) + ".csv";
        write_matrix(screen_file, screen_matrix);

                }
           }
        }

    NumericMatrix cancer_no_smk(70, 70);
    NumericMatrix cancer_smk(70, 70);
    
    for (int i = 0; i < 70; i++) {
      for (int j = 0; j < 70; j++) {
	cancer_no_smk(i, j) = mean(sim_cohort_data[i].cancer_results_no_smk(_, j));
	cancer_smk(i, j) = mean(sim_cohort_data[i].cancer_results_smk(_, j));
      }
    }

        string file_name = "Results/cancer/no_infections_" + to_string(vac_scenario) + "_" + to_string(seed) + "_" + to_string(internal_sim_N) + ".csv";
        write_matrix(file_name, no_infections);
        file_name = "Results/cancer/no_infections_thirteen_before_21_" + to_string(vac_scenario) + "_" + to_string(seed) + "_" + to_string(internal_sim_N) + ".csv";
        write_matrix(file_name, no_infections_five_before);

        file_name = "Results/cancer/mortality_denominator_" + to_string(vac_scenario) + "_" + to_string(seed) + "_" + to_string(internal_sim_N) + ".csv";
        write_matrix(file_name, mort_denom);

  file_name = "Results/cancer/cancer_smk_" + to_string(vac_scenario) + "_" + to_string(seed) + "_" + to_string(internal_sim_N) + ".csv";
        write_matrix(file_name, cancer_smk);
  file_name = "Results/cancer/cancer_no_smk_" + to_string(vac_scenario) + "_" + to_string(seed) + "_" + to_string(internal_sim_N) + ".csv";
        write_matrix(file_name, cancer_no_smk);

	return 0;
}
