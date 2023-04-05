#pragma once

#include "gbbs/gbbs.h"

namespace gbbs{

constexpr const size_t _bs_union_base = 32;
constexpr const size_t _seq_union_thresh = 2048;


template<class Seq, class F, class L, class R>
sequence<uintE> seq_union_full(const Seq& A, const Seq& B, F& f, L& l, R& r) {
	size_t nA = A.size(), nB = B.size();
	auto C = sequence<uintE>::uninitialized(nA+nB);
	size_t i = 0, j = 0;
	size_t ct = 0;
	while (i < nA && j < nB) {
		const uintE& u = A[i];
		const uintE& v = B[j];
		if (u == v) {
			if (f(u)){
				C[ct] = u;
				ct++;
			}
			i++;
			j++;
		} else if (u < v) {
			if (l(u)){
				C[ct] = u;
				ct++;
			}
			i++;
		} else {
			if (r(v)){
				C[ct] = v;
				ct++;
			}
			j++;
		}
	}
	while (i < nA){
		if (l(A[i])){
			C[ct] = A[i];
			ct++; 
		}
		i++;
	}
	while (j < nB){
		if (r(B[j])){
			C[ct] = B[j];
			ct++;
		}
		j++;
	}
	auto output = parlay::tabulate(ct, [&](size_t i){return C[i];});
	return output;
}

template<class Seq, class F, class L, class R>
sequence<uintE> seq_union(const Seq& A, const Seq& B, F& f, L& l, R& r) {
	size_t nA = A.size(), nB = B.size();
	size_t ct = 0, j = 0;
	auto C = sequence<uintE>::uninitialized(nA+nB);
	for (size_t i = 0; i < nA; i++){
		const uintE& u = A[i];
		size_t mB = parlay::binary_search(B, u, std::less<uintE>());
		if (mB < nB){
			auto filtered1 = parlay::filter_out(B.cut(j, mB), C.cut(ct, ct+mB-j), r);
			ct += filtered1;
			j = mB;
			if (u == B[mB]){
				if (f(u)){
					C[ct] = u; ct++; j++;
				} else {
					j++;
				}
			}
		} else{
			auto filtered1 = parlay::filter_out(B.cut(j, nB), C.cut(ct, ct+nB-j), r);
			ct += filtered1; j = nB;
			auto filtered2 = parlay::filter_out(A.cut(i, nA), C.cut(ct, ct+nA-i), l);
			ct += filtered2;
			break;
		}
	}
	auto filtered3 = parlay::filter_out(B.cut(j, nB), C.cut(ct, ct+nB-j), r);
	ct += filtered3; j = nB;
	auto output = parlay::tabulate(ct, [&](size_t i){return C[i];});
	return output;
}

template<class Seq, class F, class L, class R>
sequence<uintE> union_flr_par(const Seq& A, const Seq& B, 
	F& f = [](uintE x){return 1;}, L& l = [](uintE x){return 1;}, R& r = [](uintE x){return 1;}){
	size_t nA = A.size();
	size_t nB = B.size();
	size_t nR = nA + nB;
	if (nR < _seq_union_thresh) {  // handles (small, small) using linear-merge
		return seq_union_full(A, B, f, l, r);
	} else if (nB < nA) {
		return union_flr_par(B, A, f, r, l);
	} else if (nA < _bs_union_base){
		return seq_union(A, B, f, l, r);
	} else {
		size_t mA = nA / 2;
		size_t mB = parlay::binary_search(B, A[mA], std::less<uintE>());
		sequence<uintE> m_left, m_right;
		par_do(
			[&]() {m_left = union_flr_par(A.cut(0, mA), B.cut(0, mB), f, l, r);},
			[&]() {m_right = union_flr_par(A.cut(mA, nA), B.cut(mB, nB), f, l, r);}
		);
		m_left.append(m_right);
		return m_left;
	}
}

}