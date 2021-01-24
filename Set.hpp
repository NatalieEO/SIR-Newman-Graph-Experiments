/* 
 * File:   Set.hpp
 * Author: maleq
 *
 * Created on November 11, 2008, 4:35 PM
 */

#ifndef _SET_HPP
#define	_SET_HPP 

#include "utility.hpp"
				 
typedef  unsigned long int  Vertex;
typedef  unsigned long int  Vsize; 


inline int  Compare(const void *a, const void *b) 
{
	return *(Vertex *)a < *(Vertex *)b;
}


class Set {
public:
	Vsize   size;
	Vsize   maxsize;
	Vertex *set;

	Set() {
		size = 0; 
		maxsize = 0;
		set = NULL;
	}
	
	~Set() {
		if (set!=NULL) delete [] set;
	}

	Set(Vsize s) {
		size = 0; 
		maxsize = s; 
		set = new Vertex[s];
	}
	
	void  init(Vsize s) {
		if (s <= 0) return;
		maxsize = s; 
		set = new Vertex[s];
	}
	
	void decSize(){
		if(size == 0){
			cout << "Size already zero";
		}
		else
			size--;
	}
	
	void clear() {
		size = 0;
	}
	
	void destroy() {
		if (set!=NULL) delete [] set;
		set = NULL;
		size = maxsize = 0;
	}
	
	Vsize Size() const {
		return size;
	}
	
		// check if element is a member of the set -- linear search in unsorted case
	bool member(Vertex element)
	{
		return LSearch(set, size, element);
	}
	
		// Insert an element; NO boundary check - it is caller's responsibility
	void insert(Vertex element) {
		set[size] = element; 
		size++;
	}

		// Insert an element; but, dynamically increase the size of no empty space
	void Dinsert(Vertex element) { //zhao zhao: dynamic expanding set size
		if(size == maxsize){
			maxsize *= 2;
			Vertex* tmp = set;
			set = new Vertex[maxsize];
			memcpy(set, tmp, sizeof(Vertex)*size);
			delete [] tmp;		
		}
		set[size] = element; 
		size++;
	}
	
	Vertex &operator[](Vsize i) {
		return set[i];
	}
	
	int operator()(Vertex v) {		// membership function
		return BSearchD(set, 0, int(size)-1, v);  
	}
	
	void  sort() {
		qsort(set, size, sizeof(Vertex), Compare);
	}
	
	void  finalize() {
        qsort(set, size, sizeof(Vertex), Compare);
		if (size < maxsize) {
			Vertex *tmpset = set;
			set = new Vertex[size];
			for (Vsize i=0; i<size; i++) 
				set[i] = tmpset[i];
			delete [] tmpset;
		}
	}

	
	// both A and B must be sorted in descending order before calling intersect
	void intersect(Vsize i, Set &A, Set &B) {
        Vsize  k;
		size = 0;
		for (k=0; k<B.size && i<A.size; i++) {
			while (k<B.size && B.set[k] > A.set[i])  k++;
			if (k<B.size && B.set[k] == A.set[i]) insert(A.set[i]);
		}
	}
			
	// both A and B must be sorted in descending order before calling intersect
	void intersect(Set &A, Set &B) {
        Vsize  i, k;
		size = 0;
		
		for (i=0, k=0; k<B.size && i<A.size; i++) {
			while (k<B.size && B.set[k] > A.set[i])  k++;
			if (k<B.size && B.set[k] == A.set[i]) insert(A.set[i]);
		}
	}
	
	// both A and B must be sorted in descending order before calling intersect
	void intersect(Set &B) {
        Vsize  i, k;
		Vsize tsize = size;
		size = 0;
		
		for (i=0, k=0; k<B.Size() && i<tsize; i++) {
			while (k<B.Size() && B[k] > set[i])  k++;
			if (k<B.Size() && B[k] == set[i]) insert(set[i]);
		}
	}

	// both A and B must be sorted in descending order before calling set_union
	void set_union(Set &A, Set &B) {
        Vsize  i, k;
		size = 0;
		
		for (i=0, k=0; k<B.Size() && i<A.Size(); i++) {
			while (k<B.Size() && B[k] > A[i])  insert(B[k++]);
			if (k>=B.Size() || B[k] != A[i])  insert(A[i]);
		}
		while (i<A.Size())  insert(A[i++]);
		while (k<B.Size())  insert(B[k++]);
	}
	
	// both A and B must be sorted in descending order before calling set_union
	void set_union(Set &B) {
        Vsize  i, k;
		Vertex *tset = set;
		Vsize tsize = size; 
		
		size = 0;
		init(tsize+B.Size());
		
		for (i=0, k=0; k<B.Size() && i<tsize; i++) {
			while (k<B.Size() && B[k] > tset[i])  insert(B[k++]);
			if (k>=B.Size() || B[k] != tset[i])  insert(tset[i]);
		}
		while (i<tsize)  insert(tset[i++]);
		while (k<B.Size())  insert(B[k++]);
		
		if (tset) delete [] tset;
	}
		
	void print() {
		cout<< "Size "<< size <<": "; 
		for (int i=0; i<size; i++) 
			cout << set[i] << " "; 
		cout<<endl;
	}

	void intersect_idx(Set &A, Set &B) {
        Vsize  i, k;
		size = 0;
		
		for (i=0, k=0; k<B.size && i<A.size; i++) {
			while (k<B.size && B[k] > A[i])  k++;
			if (k<B.size && B[k] == A[i]) insert(i);
		}
	}
};

#endif	/* _SET_HPP */

