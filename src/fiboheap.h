#include <iostream>

using namespace std;

template <class V> class FiboHeap;

template <class V> class FiboNode {
private:
	FiboNode<V>* prev;
	FiboNode<V>* next;
	FiboNode<V>* child;
	FiboNode<V>* parent;
	V value;
	int degree;
	bool marked;
public:
	friend class FiboHeap<V>;
	FiboNode<V>* getPrev() {return prev;}
	FiboNode<V>* getNext() {return next;}
	FiboNode<V>* getChild() {return child;}
	FiboNode<V>* getParent() {return parent;}
	V getValue() {return value;}
	bool isMarked() {return marked;}
	bool hasChildren() {return child;}
	bool hasParent() {return parent;}
};

template <class V> class FiboHeap {
protected:
	// circular doubly linked list
	FiboNode<V>* heap;
	static const int MAXDEGREE = 5000;
public:
	FiboHeap() {
		heap = _empty();
	}
	~FiboHeap() {
		if(heap)
			_deleteAll(heap);
	}
	FiboNode<V>* insert(V value) {
		FiboNode<V>* ret = _singleton(value);
		heap = _merge(heap,ret);
		return ret;
	}
	void merge(FiboHeap& other) {
		heap = _merge(heap,other.heap);
		other.heap = _empty();
	}

	bool isEmpty() {
		return heap == NULL;
	}

	V getMinimum() {
		return heap->value;
	}

	V removeMinimum() {
		FiboNode<V>* old = heap;
		heap = _removeMinimum(heap);
		V ret = old->value;
		delete old;
		return ret;
	}

	void decreaseKey(FiboNode<V>* n,V value) {
		heap = _decreaseKey(heap,n,value);
	}

	void update(FiboNode<V>* n) {
		if (n == heap)
			heap = _getMinIn(heap);
		heap = _update(heap,n);
	}

	V getChild() {
		return heap->child->getValue();
	} 

private:
	FiboNode<V>* _empty() {
		heap = NULL;
		return heap;
	}

	FiboNode<V>* _singleton(V value) {
		FiboNode<V>* n = new FiboNode<V>;
		n->value = value;
		n->prev = n;
		n->next = n;
		n->degree = 0;
		n->marked = false;
		n->child = NULL;
		n->parent = NULL;
		return n;
	}

	// a and b are circular doubly linked lists
	FiboNode<V>* _merge(FiboNode<V>* a,FiboNode<V>* b) {
		if(a == NULL)	return b;
		if(b == NULL)	return a;
		FiboNode<V>* an = a->next;
		FiboNode<V>* bp = b->prev;
		a->next = b;
		b->prev = a;
		an->prev = bp;
		bp->next = an;
		return ((*(a->value) < *(b->value)) ? a : b);
	}

	void _deleteAll(FiboNode<V>* n) {
		if(n == NULL)	return;
		FiboNode<V>* p = n;
		do {
			FiboNode<V>* d = p;
			p = p->next;
			_deleteAll(d->child);
			delete d;
		} while(p != NULL && p != n);
	}

	void _addChild(FiboNode<V>* parent,FiboNode<V>* child) {
		child->prev = child;
		child->next = child;
		child->parent = parent;
		parent->degree++;
		parent->child = _merge(parent->child,child);
	}

	void _unMarkAndUnParentAll(FiboNode<V>* n) {
		if(n == NULL)	return;
		FiboNode<V>* p = n;
		do {
			p->marked = false;
			p->parent = NULL;
			p = p->next;
		}while(p != n);
	}

	FiboNode<V>* _getMinIn(FiboNode<V>* n) {
		FiboNode<V>* min = n;
		FiboNode<V>* p = n;
		do {
			if(*(p->value) < *(min->value))
				min = p;
			p = p->next;
		} while(p != n);
		return min;
	}

	FiboNode<V>* _removeMinimum(FiboNode<V>* n) {
		// meld
		_unMarkAndUnParentAll(n->child);
		if(n->next == n) {
			n = n->child;
		} else {
			n->next->prev = n->prev;
			n->prev->next = n->next;
			n = _merge(n->next,n->child);
		}
		if(n == NULL)	return n;
		
		// consolidate
		FiboNode<V>* trees[MAXDEGREE] = {NULL};
		while(true) {
			if(trees[n->degree] != NULL) {
				FiboNode<V>* t = trees[n->degree];
				if(t == n)	break;
				trees[n->degree] = NULL;
				if(*(n->value) < *(t->value)) {
					t->prev->next = t->next;
					t->next->prev = t->prev;
					_addChild(n,t);
				} else {
					t->prev->next = t->next;
					t->next->prev = t->prev;
					if(n->next == n) {
						t->next = t->prev = t;
						_addChild(t,n);
						n = t;
					} else {
						n->prev->next = t;
						n->next->prev = t;
						t->next = n->next;
						t->prev = n->prev;
						_addChild(t,n);
						n = t;
					}
				}
				continue;
			} else {
				trees[n->degree] = n;
			}
			n = n->next;
		}
		// update min
		return _getMinIn(n);
	}

	FiboNode<V>* _cut(FiboNode<V>* heap,FiboNode<V>* n) {
		if(n->next == n) {
			n->parent->child = NULL;
		} else {
			n->next->prev = n->prev;
			n->prev->next = n->next;
			if (n->parent->child == n)
				n->parent->child = _getMinIn(n->next);
		}
		n->next = n;
		n->prev = n;
		n->marked = false;
		n->parent->degree--;
		n->parent = NULL;
		return _merge(heap,n);
	}

	FiboNode<V>* _decreaseKey(FiboNode<V>* heap, FiboNode<V>* n, V value) {
		if(*(n->value) < *value){
			return heap;
		}
		n->value = value;
		// update structure
		if (n->parent == NULL && *(n->value) < *(heap->value))
			heap = n;
		if(n->parent != NULL && *(n->value) < *((n->parent)->value)) {
			FiboNode<V>* parent = n->parent;
			heap = _cut(heap,n);
			while(parent->parent != NULL && parent->marked) {
				FiboNode<V>* tmp = parent->parent;
				heap = _cut(heap, parent);
				parent = tmp;
			}
			if(parent != NULL && parent->parent != NULL)
				parent->marked = true;
		}
		return heap;
	}
	
	// only appliable for increase key
	FiboNode<V>* _update(FiboNode<V>* heap,FiboNode<V>* n) {
		if (n->child == NULL)	return heap;
		FiboNode<V>* kids = n->child;
		kids = _getMinIn(kids);
		if (*(n->value) < *(kids->value))
			return heap;
		n->child = NULL;
		return _merge(heap, kids);
	}

};
