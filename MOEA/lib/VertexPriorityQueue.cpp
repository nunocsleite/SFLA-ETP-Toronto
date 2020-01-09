
#include "VertexPriorityQueue.h"


VertexPriorityQueue::VertexPriorityQueue(int nvertices) : handles(nvertices+1) { }

void VertexPriorityQueue::push(int vertex, int priority) {
	heap_data data(vertex, priority);
	handle_t handle = pq.push(data);
    // store handle
	handles[vertex] = handle;
}

VertexPriorityQueue::heap_data const& VertexPriorityQueue::top() const {
	return pq.top();
}

void VertexPriorityQueue::pop() {
	pq.pop();
}

void VertexPriorityQueue::update(int vertex, int newPriority) {
	// get handle
	handle_t handle = handles[vertex];
	(*handle).priority = newPriority;
	pq.update(handle);
}

bool VertexPriorityQueue::empty() const {
	return pq.empty();
}

void VertexPriorityQueue::printQueue() const {
		cout << "Priority Queue: ordered iteration" << endl;
	for (PQ::ordered_iterator it = pq.ordered_begin(); it != pq.ordered_end(); ++it)
		cout << *it << " "; 
	cout << endl;
}



//
//int main() {
//
//	int N = 10;
//
//	VertexPriorityQueue pq(N);
//
//
//	cin.get();
//
//	return 0;
//}
