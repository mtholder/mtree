#if !defined(__TREE_H__)
#define __TREE_H__
#include <cassert>
#include <stack>
#include <climits>
#include <iostream>
#include <vector>
#include <algorithm>
#include "mt_log.h"
#include "mt_model_description.h"
namespace mt {
class PartitionedMatrix;
class LeafCharacterVector;
class LeafWork;
class InternalNodeWork;

class Node {
    public:
        Node()
            :parent(nullptr),
             leftChild(nullptr),
             rightSib(nullptr),
             number(UINT_MAX),
             scoreFlag(true),
             edgeLen(-1.0) {
        }
        std::vector<Node *> GetChildren() const {
            std::vector<Node *> c;
            Node * curr = leftChild;
            while (curr) {
                c.push_back(curr);
                curr = curr->rightSib;
            }
            return c;
        }
        void SetNumber(unsigned i) {
            this->number = i;
        }
        unsigned GetNumber() const {
            return this->number;
        }
        void SetFlag(bool v) {
          this->scoreFlag = v;
        }
        bool GetFlag() {
          return this->scoreFlag;
        }
        void AddChild(Node *c, double newEdgeLen) {
            assert(c);
            c->parent = this;
            Node *r = this->GetLastChild();
            if (r == nullptr) {
                this->leftChild = c;
            } else {
                r->rightSib = c;
            }
            assert(c->rightSib == nullptr);
            c->rightSib = nullptr;
            c->SetEdgeLen(newEdgeLen);
        }
        Node * GetLastChild() {
            if (this->leftChild == nullptr) {
                return nullptr;
            }
            Node * c = this->leftChild;
            while (c->rightSib) {
                c = c->rightSib;
            }
            return c;
        }
        Node * GetRoot() {
          Node * r = this;
          while(r) {
            if(r->parent == nullptr) return r;
            r = r->parent;
          }
        }
        void SetEdgeLen(double e) {
            this->edgeLen = e;
        }
        double GetEdgeLen() const {
            return this->edgeLen;
        }
        void SetVEdgeLen(double v) {
            this->vEdgeLen = v;
        }
        double GetVEdgeLen() const {
            return this->vEdgeLen;
        }
        bool IsLeaf() const {
            return this->leftChild == nullptr;
        }
        bool IsLeftChild() const {
            return this->rightSib != nullptr;
        }
        void SetData(unsigned i, void * d) {
            while (this->data.size() <= i) {
                this->data.push_back(nullptr);
            }
            this->data[i] = d;
        }
        void SetWork(unsigned i, void * d) {
            while (this->work.size() <= i) {
                this->work.push_back(nullptr);
            }
            this->work[i] = d;
        }
        void * GetWork(unsigned i) {
            return work[i];
        }
        void * GetData(unsigned i) {
            return data[i];
        }
        void write(std::ostream & out) const;
        Node * parent;
        Node * leftChild;
        Node * rightSib;
    private:
        unsigned number;
        double edgeLen;
        double vEdgeLen;
        bool scoreFlag;
        std::vector<void *> data;
        std::vector<void *> work;
        friend class Arc;
};

class Tree {
    public:
        std::size_t GetNumLeaves() const {
            return leaves.size();
        }
        std::size_t GetNumNodes() const {
            return nodes.size();
        }
        void SetRoot(Node *r) {
            assert(r);
            assert(r->GetNumber() >= this->GetNumLeaves());
            this->root = r;
        }
        Node * GetRoot() {
            return this->root;
        }
        const Node * GetRoot() const {
            return this->root;
        }
        Node * GetNode(unsigned i) {
            return &(this->nodes.at(i));
        }
        Node * GetLeaf(unsigned i) {
            return this->leaves[i];
        }
        const Node * GetLeaf(unsigned i) const {
            return this->leaves[i];
        }
        // tree method to check if node is a leaf
        bool isLeaf2(Node * q) {
            return (q->GetNumber() < leaves.size());
        }
        Tree(unsigned numNodes, unsigned numLeaves)
            :nodes(numNodes),
             root(nullptr) {
            initPointers(numLeaves);
        }
        Tree(const Tree &other)
            :nodes(other.nodes.size()),
            root(nullptr) {
            initPointers(other.leaves.size());
            copyTopology(other);
        }
        void copyTopology(const Tree & other);
        bool isNodeConnected(Node *n, int id);
        void write(std::ostream & out) const {
            root->write(out);
            out << ";\n";
        }
        void TreeDebug() const;

    private:
        void initPointers(unsigned numLeaves) {
          const unsigned numNodes = nodes.size();
          assert(numLeaves < numNodes);
          leaves.resize(numLeaves);
          for (auto i = 0U; i < numNodes; ++i) {
              Node & node = this->nodes[i];
              node.SetNumber(i);
              if (i < numLeaves) {
                  leaves[i] = &node;
              }
          }
        }
        std::vector<Node> nodes;
        Node * root;
        std::vector<Node *> leaves;
};

class CharModel;

class InternalNodeWork {
    public:
        InternalNodeWork(std::size_t numChars, unsigned numStates, unsigned numRates)
            :claAtNdFromNd(numStates*numRates*numChars),
             claAtNdFromPar(numStates*numRates*numChars),
             claAtParFromNd(numStates*numRates*numChars),
             claAtParFromPar(numStates*numRates*numChars),
             nChars(numChars) {
            //_DEBUG_VAL(numChars);
            //_DEBUG_VAL(numStates);
            //_DEBUG_VAL(numRates);
            //_DEBUG_VEC(claAtNdFromNd);
            //_DEBUG_VEC(claAtNdFromPar);
            //_DEBUG_VEC(claAtParFromNd);
            //_DEBUG_VEC(claAtParFromPar);
        }
        std::size_t GetLenCLA() const {
            return claAtNdFromNd.size();
        }
    std::vector<double> claAtNdFromNd;
    std::vector<double> claAtNdFromPar;
    std::vector<double> claAtParFromNd;
    std::vector<double> claAtParFromPar;
    const std::size_t nChars;
};

class LeafWork: public InternalNodeWork {
    public:
        LeafWork(std::size_t numChars, unsigned numStateCodes, unsigned numStates, unsigned numRates)
            :InternalNodeWork(numChars, numStates, numRates),
            summed(numStateCodes*numStates*numRates){
                //_DEBUG_VAL(numStateCodes);
                //_DEBUG_VEC(summed);
            }
        double * GetCLAElements() {
            return &(summed[0]);
        }
    std::vector<double> summed;
};


class Arc {
    public:
        Arc(Node * fromNd, Node * toNd)
            :fromNode(fromNd),
            toNode(toNd),
            edgeLenPtr(nullptr),
            fromIsChild(false) {
                if (toNd) {
                    if (toNode->parent == fromNode) {
                        fromIsChild = false;
                        edgeLenPtr = &(toNode->edgeLen);
                    } else {
                        assert(fromNode->parent == toNode);
                        fromIsChild = true;
                        edgeLenPtr = &(fromNode->edgeLen);
                    }
                } else { /* not really a "child" */
                    fromIsChild = true;
                }
            }
        double GetEdgeLen() const {
            return *edgeLenPtr;
        }
        void SetEdgeLen(double x) {
            *edgeLenPtr = x;
        }
        bool IsFromLeaf() const {
            return fromNode->IsLeaf();
        }
        bool IsToLeaf() const {
            return toNode->IsLeaf();
        }
        const LeafCharacterVector * GetFromNdData(unsigned partIndex) {
            assert(IsFromLeaf());
            return static_cast<const LeafCharacterVector *>(fromNode->GetData(partIndex));
        }
        LeafWork * GetFromNdLeafWork(unsigned partIndex) {
            assert(IsFromLeaf());
            return static_cast<LeafWork *>(fromNode->GetWork(partIndex));
        }
        InternalNodeWork * GetFromNdIntWork(unsigned partIndex) {
            return static_cast<InternalNodeWork *>(fromNode->GetWork(partIndex));
        }
        InternalNodeWork * GetToNdIntWork(unsigned partIndex) {
            return static_cast<InternalNodeWork *>(toNode->GetWork(partIndex));
        }
        double * GetFromNdCLA(unsigned partIndex, bool crossedEdge);
        double * GetToNdCLA(unsigned partIndex, bool crossedEdge);
        std::vector<const double *> GetPrevCLAs(unsigned partIndex);
        std::size_t GetLenCLA(unsigned partIndex) {
            return GetFromNdIntWork(partIndex)->GetLenCLA();
        }
        std::size_t GetNumChars(unsigned partIndex) {
            return GetFromNdIntWork(partIndex)->nChars;
        }
        Node * fromNode;
        Node * toNode;
    private:
        double * edgeLenPtr;
        bool fromIsChild;
};

inline double * Arc::GetFromNdCLA(unsigned partIndex, bool crossedEdge) {
    if (fromIsChild) {
        InternalNodeWork * work = GetFromNdIntWork(partIndex);
        if (crossedEdge) {
            return &(work->claAtParFromNd[0]);
        } else {
            return &(work->claAtNdFromNd[0]);
        }
    } else {
        InternalNodeWork * work = GetToNdIntWork(partIndex);
        if (crossedEdge) {
            return &(work->claAtNdFromPar[0]);
        } else {
            return &(work->claAtParFromPar[0]);
        }
    }
}

inline double * Arc::GetToNdCLA(unsigned partIndex, bool crossedEdge) {
    if (fromIsChild) {
        InternalNodeWork * work = GetFromNdIntWork(partIndex);
        if (crossedEdge) {
            return &(work->claAtNdFromPar[0]);
        } else {
            return &(work->claAtParFromPar[0]);
        }
    } else {
        InternalNodeWork * work = GetToNdIntWork(partIndex);
        if (crossedEdge) {
            return &(work->claAtParFromNd[0]);
        } else {
            return &(work->claAtNdFromNd[0]);
        }
    }
}

inline std::vector<const double *> GetSurroundingCLA(Node * fromNode, Node * avoid, unsigned partIndex) {
    std::vector<const double *> pcla;
    std::vector<Node *> c = fromNode->GetChildren();
    for (auto i : c) {
        if (i != avoid) {
            InternalNodeWork * iw = static_cast<InternalNodeWork *>(i->GetWork(partIndex));
            double * ic = &(iw->claAtParFromNd[0]);
            pcla.push_back(const_cast<const double *>(ic));
        }
    }
    if (fromNode->parent && fromNode->parent != avoid) {
        InternalNodeWork * iw =  static_cast<InternalNodeWork *>(fromNode->GetWork(partIndex));
        double * ic = &(iw->claAtNdFromPar[0]);
        pcla.push_back(const_cast<const double *>(ic));
    }
    return pcla;
}

inline std::vector<const double *> Arc::GetPrevCLAs(unsigned partIndex) {
    std::vector<const double *> pcla;
    if (fromIsChild) {
        assert(fromNode->leftChild);
        std::vector<Node *> c = fromNode->GetChildren();
        for (auto i : c) {
            InternalNodeWork * iw =  static_cast<InternalNodeWork *>(i->GetWork(partIndex));
            double * ic = &(iw->claAtParFromNd[0]);
            pcla.push_back(const_cast<const double *>(ic));
        }
    } else {
        return GetSurroundingCLA(fromNode, toNode, partIndex);
    }
    return pcla;
}


} //namespace
#endif
