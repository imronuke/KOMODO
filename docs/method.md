# Theory and Background

Since ADPRES uses the current available methods, rather than re-explaining the methods, we would like to just mention the references of the methods implemented in the ADPRES.

The ADPRES development began in 2017 and was motivated by many hurdles to obtain such similar codes in the region where the author was working. The initial version of ADPRES employed Nodal Expansion Method (NEM) based on the response matrix formulations [1,2,3,4]. The transverse integrated leakages were approximated by the Quadratic Transverse Leakage Approximation (QTLA)[1]. The transient and and thermal module were also added to enable ADPRES solving time-dependent problems with thermal-hydraulics (T-H) feedback. The transient diffusion equation was solved with fully implicit method, and T-H solutions were obtained by solving mass and energy conservation equations in an enclosed channel following NODAL 3 computer code [5] with slight modifications. This version of ADPRES was published in Annals of Nuclear Energy in 2019 [6]. However, this version encountered slow performance compared to other modern nodal simulators.

Thus, in the early 2020, the ADPRES was revamped to implement CMFD acceleration with two-node problems and non-linear iteration procedures [7] for the sake of rapid calculations. By implementing CMFD acceleration, the current ADPRES version is able to solve time-dependent problems 10-20 times faster than previous version of ADPRES. Initially, ADPRES implemented Polynomial Nodal Method [8,5], then upgraded to Semi-Analytic Nodal Method [9,10] for better accuracy and the neutron precursor equations are solved analytically [11]. Also, in the transient calculations, the delayed terms and other terms that do not appear in static calculations are included with transverse leakages in the calculation of transverse moments to further save memory storage [12]. Users also have option to perform exponential flux [13] transformation and to set theta value for time-dependent problems. It is also possible to perform calculations with branched cross sections (i.e. set of cross sections for various TH parameters, boron concentration, and control rod position). The branched cross section library format was adopted from [14] and had been tested UO-MOX PWR transient benchmark [14].

# Acknowledgement

This ADPRES development would have been impossible without the God's Mercy and the works done by those mentioned in the references. We would like to thank them and other people who contributed on their works. We also would like to thank to other people who directly or indirectly contributed to this work:

* Dr. Ali Al Naqbi, Abu Dhabi Polytechnic
* Dr. Anthony Hechanova, Abu Dhabi Polytechnic
* Prof. Nam Zin Cho, KAIST
* Dr. Alexander Agung, Gadjah Mada University
* Dr. Andang Widiharto, Gadjah Mada University
* Liem Peng Hong, PhD, NAIS and Tokyo City University
* Donny Hartanto, Phd
* GNU Fortran and Intel Fortran developer teams.
* All my colleagues and friends.

# References

1. Finnemann, H., Bennewitz F. and Wagner M. R., (1977) Interface current techniques for multidimensional reactor calculations, Atomkernenergie, Vol. 30, pp. 123-128.

2. Bennewitz F., Finnemann H. and Moldaschl H.(1975)  Solution of the multidimensional neutron diffusion equa-tion by nodal expansion. Proc. Conf. Computational Methods in Nuclear Engineering, p. 1-99, CONF-750413.

3. Lawrence, R.D., (1986) Progress in nodal methods for the solution of the neutron diffusion and transport equations, Progress in Nuclear Energy, Vol. 17, No.3, pp. 271-301.

4. Okumura, K., (1998) MOSRA-Light: High speed three-dimensional nodal diffusion code for vector computers, JAEA-Research 98-025. (in Japanese)

5. Liem P.H., et al., (2010) NODAL3: A Nodal Neutron Diffusion Code Version 2.0 User Guides (unpublihsed)

6. Imron, M. (2019). Development and verification of open reactor simulator ADPRES. Annals of Nuclear Energy, 133, 580–588. https://doi.org/10.1016/j.anucene.2019.06.049

7. Smith, K. S. (1984) Nodal Method Storage Reduction by Nonlinear Iteration. Transactions of American Nuclear Society 44, 265.

8. Zimin V.G. and Ninokata, H., (1997) Nonlinear Iteration Procedure Based on Legendre Polynomials, Trans. Am. Nucl. Soc., 76, 162.

9. Zimin, V. G., & Ninokata, H. (1998). Nodal neutron kinetics model based on nonlinear iteration procedure for LWR analysis. Annals of Nuclear Energy, 25(8), 507–528. https://doi.org/10.1016/S0306-4549(97)00078-9

10. Zimin V.G. and Ninokata, H., (1997) Polynomials ans Semi-Analytic Nodal Methods For Non-Linear Iteration Procedures, PYHSOR-98, 2, 994.

11. Stacey, W. M. (1969) Space-Time Nuclear Reactor Kinetics. Academic Press, New York.

12. Engrand, P. R., Maldonado, G. I., A1-Chalabi, R. M. and Turinsky, P. J. (1992) Non-Linear Iterative Strategy for NEM Refinement and Extension. Transactions of American Nuclear Society 65, 221.

13. Hendricks, J.S., “Finite difference solution of the time dependent neutron group diffusion equations”, Thesis, Department of Nuclear Engineering, Massachusetts Institute of Technology, MITNE-176 (1975).

14. Kozlowski, T., Downar, T., 2007. PWR MOX/UO2 Core Transient Benchmark Final Report. Retrieved from: https://www.oecd-nea.org/science/wprs/MOX-UOX-transients/benchmark_documents/final_report/mox-uo2-bench.pdf
