
# include <SOFT/TROTTER_prop.hpp>
# include <SOFT/constants.hpp>
# include <SOFT/fileops.hpp>
# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <numeric>
# include <fftw3.h>





int main()	{
	const int wp_size = std::accumulate(std::begin(
				DYN_CONST::grid_pt_vect),std::end(DYN_CONST::grid_pt_vect)
			,1,std::multiplies<int>());
	std::vector<fftw_complex> psi_x(wp_size);
	std::vector<int> line_extract_no_vect(wp_size);
	for (int i=0; i<wp_size; i++)	{
		line_extract_no_vect[i] = i+1;
	}
	/**------------------------------------**/
	/** Extracting data from potential and initial wave packet files **/
	std::vector<double> pot_vect;
	std::vector<double> psi_x_init_real;
	if (DYN_CONST::initfiletype == "dat")	{
		pot_vect = get_all_sp_lines(line_extract_no_vect,DYN_CONST::potfile);
		psi_x_init_real = get_all_sp_lines(line_extract_no_vect,DYN_CONST::initfile); 
		for (size_t i=0; i<psi_x.size(); i++)	{
			psi_x [i][0] = psi_x_init_real[i];
			psi_x [i][1] = 0.;
		}
	}
	else if (DYN_CONST::initfiletype == "bin")	{
		std::vector<double> dummyvect;
		pot_vect = rtrn_vec_from_bin(dummyvect,DYN_CONST::potfile,1,wp_size);
		psi_x_init_real = rtrn_vec_from_bin(dummyvect,DYN_CONST::initfile,1,wp_size);
		for (size_t i=0; i<psi_x.size(); i++)	{
			psi_x [i][0] = psi_x_init_real[i];
			psi_x [i][1] = 0.;
		}
	}
	/**------------------------------------**/
	std::vector<std::vector<double>> p_mult_grid_vect; 
	/*----------------------------*/
	std::vector<double> delta_vect(DYN_CONST::grid_pt_vect.size());
	for (size_t i=0; i<delta_vect.size(); i++)	{
		delta_vect[i] = (DYN_CONST::grid_up_lim_vect[i]-
				DYN_CONST::grid_low_lim_vect[i])/(DYN_CONST::grid_pt_vect[i]-1);
	}
	/*----------------------------*/
	for (size_t i=0; i<delta_vect.size(); i++)	{
		std::vector<double> p_grid_1d(DYN_CONST::grid_pt_vect[i]);
		onedim_momentum_grid(p_grid_1d,DYN_CONST::grid_pt_vect[i],delta_vect[i]);
		p_mult_grid_vect.push_back(p_grid_1d);
		empty_swap(p_grid_1d);
	}
	/*----------------------------*/
	std::vector<fftw_complex> trotter_pot = get_pot_trotter_part(pot_vect);
	std::vector<fftw_complex> trotter_ke = get_ke_trotter_part(p_mult_grid_vect,DYN_CONST::freq_vect);
	int sizevect[INIT_STATE_CONST::NDIM_4_DYNAMICS];
	int rank = INIT_STATE_CONST::NDIM_4_DYNAMICS; 
	for (int i=0; i<INIT_STATE_CONST::NDIM_4_DYNAMICS; i++)	{
		sizevect[i] = DYN_CONST::grid_pt_vect[i];
	}
	fftw_plan plan_forward;
	fftw_plan plan_backward;
	plan_Multidim_ForwardFFT(plan_forward,psi_x.size(),sizevect,rank,"in place");
	plan_Multidim_BackwardFFT(plan_backward,psi_x.size(),sizevect,rank,"in place");
	int counter = 0;
	std::vector<double> psi_real_imag_abs_vect(3*wp_size);		// vector containing real,imaginary and absolute  values of psi
	rtrn_real_imag_norm(psi_real_imag_abs_vect,psi_x);
	save_binary(psi_real_imag_abs_vect,DYN_CONST::save_outbinfile,"append");
	save_binary(psi_real_imag_abs_vect,DYN_CONST::restart_outbinfile,"append");
	psi_real_imag_abs_vect.assign(psi_real_imag_abs_vect.size(),0);
	for (int i=0; i<DYN_CONST::max_step+1; i++)	{
		generate_fftw_complex_product(trotter_pot,psi_x);
		execute_inplaceFFT(psi_x,"FORWARD",plan_forward);		
		generate_fftw_complex_product(trotter_ke,psi_x);
		execute_inplaceFFT(psi_x,"BACKWARD",plan_backward);
		generate_fftw_complex_product(trotter_pot,psi_x);
		counter ++;
		if (counter% DYN_CONST::step_4_file == 0)	{
			rtrn_real_imag_norm(psi_real_imag_abs_vect,psi_x);
			save_binary(psi_real_imag_abs_vect,DYN_CONST::save_outbinfile,"append");
			psi_real_imag_abs_vect.assign(psi_real_imag_abs_vect.size(),0);
		}
		else if (counter% DYN_CONST::step_4_restart == 0)	{
			rtrn_real_imag_norm(psi_real_imag_abs_vect,psi_x);
			save_binary(psi_real_imag_abs_vect,DYN_CONST::restart_outbinfile,"append");
			psi_real_imag_abs_vect.assign(psi_real_imag_abs_vect.size(),0);
		}
	}
	return  EXIT_SUCCESS;
}
