#if !defined(__TREE_H__)
#define __TREE_H__
#include <cassert>
#include <stack>
#include <climits>
#include <iostream>
#include <vector>
#include "mt_log.h"
namespace mt {
class PartitionedMatrix;
class LeafCharacterVector;
class ModelDescription {
    public:
        enum AscBiasMode {
            NO_ASC_BIAS = 0,
            VAR_ONLY_NO_MISSING_ASC_BIAS = 1,
            VAR_ONLY_MISSING_ASC_BIAS = 2,
            PARS_ONLY_NO_MISSING_ASC_BIAS = 3,
            PARS_ONLY_MISSING_ASC_BIAS = 4
        };
        ModelDescription(AscBiasMode m)
            :ascBiasMode(m) {
        }
        const AscBiasMode & GetAscBiasMode() const {
            return this->ascBiasMode;
        }
    private:
        AscBiasMode ascBiasMode;
};
class LeafWork;
class InternalNodeWork;

class Node {
    public:
        Node()
            :parent(nullptr),
             leftChild(nullptr),
             rightSib(nullptr),
             number(UINT_MAX),
             edgeLen(-1.0) {
        }
        std::vector<const Node *> GetChildren() const {
            std::vector<const Node *> c;
            const Node * curr = leftChild;
            while (curr) {
                c.push_back(curr);
                curr = curr->rightSib;
            }
            return c;
        }
        std::vector<Node *> GetChildren() {
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
        void AddChild(Node *c, double edgeLen) {
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
            c->SetEdgeLen(edgeLen);
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
        void SetEdgeLen(double e) {
            this->edgeLen = e;
        }
        double GetEdgeLen() {
            return this->edgeLen;
        }
        bool IsLeaf() const {
            return this->leftChild == nullptr;
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
        void * GetWork(unsigned i) const {
            return work[i];
        }
        void * GetData(unsigned i) {
            return data[i];
        }
    private:
    public:
        Node * parent;
        Node * leftChild;
        Node * rightSib;
        unsigned number;
        double edgeLen;
        std::vector<void *> data;
        std::vector<void *> work;
};
class Tree {
    public:
        unsigned GetNumLeaves() const {
            return leaves.size();
        }
        void SetRoot(Node *r) {
            assert(r);
            assert(r->number >= this->GetNumLeaves());
            this->root = r;
        }
        Node * GetRoot() {
            return this->root;
        }
        const Node * GetRoot() const {
            return this->root;
        }
        Node * GetNode(unsigned i) {
            return &(this->nodes[i]);
        }
        Node * GetLeaf(unsigned i) {
            return this->leaves[i];
        }
        const Node * GetLeaf(unsigned i) const {
            return this->leaves[i];
        }
        Tree(unsigned numNodes, unsigned numLeaves)
            :nodes(numNodes),
             root(nullptr) {
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
    private:
        std::vector<Node> nodes;
        Node * root;
        std::vector<Node *> leaves;
};

class CharModel;

class InternalNodeWork {
    public:
        InternalNodeWork(unsigned numChars, unsigned numStates, unsigned numRates)
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
        unsigned GetLenCLA() const {
            return claAtNdFromNd.size();
        }
    std::vector<double> claAtNdFromNd;
    std::vector<double> claAtNdFromPar;
    std::vector<double> claAtParFromNd;
    std::vector<double> claAtParFromPar;
    const unsigned nChars;
};

class LeafWork: public InternalNodeWork {
    public:
        LeafWork(unsigned numChars, unsigned numStateCodes, unsigned numStates, unsigned numRates)
            :InternalNodeWork(numChars, numStates, numRates),
            summed(numStateCodes*numStates*numRates){
                _DEBUG_VAL(numStateCodes);
                _DEBUG_VEC(summed);
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
        const LeafCharacterVector * GetFromNdData(unsigned partIndex) const {
            assert(IsFromLeaf());
            return (const LeafCharacterVector *)fromNode->GetData(partIndex);
        }
        LeafWork * GetFromNdLeafWork(unsigned partIndex) const {
            assert(IsFromLeaf());
            return (LeafWork *)fromNode->GetWork(partIndex);
        }
        InternalNodeWork * GetFromNdIntWork(unsigned partIndex) const {
            return (InternalNodeWork *)fromNode->GetWork(partIndex);
        }
        InternalNodeWork * GetToNdIntWork(unsigned partIndex) const {
            return (InternalNodeWork *)toNode->GetWork(partIndex);
        }

        double * GetFromNdCLA(unsigned partIndex, bool crossedEdge) const;
        double * GetToNdCLA(unsigned partIndex, bool crossedEdge) const;
        std::vector<const double *> GetPrevCLAs(unsigned partIndex) const;
        unsigned GetLenCLA(unsigned partIndex) const {
            return GetFromNdIntWork(partIndex)->GetLenCLA();
        }
        unsigned GetNumChars(unsigned partIndex) const {
            return GetFromNdIntWork(partIndex)->nChars;
        }
        Node * fromNode;
        Node * toNode;
    private:
        double * edgeLenPtr;
        bool fromIsChild;
};

class ConstArc {
    Arc arc;
    public:
        ConstArc(const Arc & a)
            :arc(a) {
        }
        ConstArc(const Node * fromNd, const Node * toNd)
            :arc(const_cast<Node *>(fromNd), const_cast<Node *>(toNd)) {
        }
        const Node * toNode() const {
            return arc.toNode;
        }
        const Node * fromNode() const {
            return arc.toNode;
        }
        double GetEdgeLen() const {
            return arc.GetEdgeLen();
        }
        bool IsFromLeaf() const {
            return arc.IsFromLeaf();
        }
        bool IsToLeaf() const {
            return arc.IsToLeaf();
        }
        const LeafCharacterVector * GetFromNdData(unsigned partIndex) const {
            return arc.GetFromNdData(partIndex);
        }
        LeafWork * GetFromNdLeafWork(unsigned partIndex) const {
           return arc.GetFromNdLeafWork(partIndex);
        }
        InternalNodeWork * GetFromNdIntWork(unsigned partIndex) const {
            return arc.GetFromNdIntWork(partIndex);
        }
        InternalNodeWork * GetToNdIntWork(unsigned partIndex) const {
            return arc.GetToNdIntWork(partIndex);
        }
        double * GetFromNdCLA(unsigned partIndex, bool crossedEdge) const {
            return arc.GetFromNdCLA(partIndex, crossedEdge);
        }
        double * GetToNdCLA(unsigned partIndex, bool crossedEdge) const {
            return arc.GetToNdCLA(partIndex, crossedEdge);
        }
        std::vector<const double *> GetPrevCLAs(unsigned partIndex) const {
            return arc.GetPrevCLAs(partIndex);
        }
        unsigned GetLenCLA(unsigned partIndex) const {
            return arc.GetLenCLA(partIndex);
        }
        unsigned GetNumChars(unsigned partIndex) const {
            return arc.GetNumChars(partIndex);
        }
};

inline double * Arc::GetFromNdCLA(unsigned partIndex, bool crossedEdge) const {
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

inline double * Arc::GetToNdCLA(unsigned partIndex, bool crossedEdge) const {
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

inline std::vector<const double *> GetSurroundingCLA(const Node * fromNode,
                                                     const Node * avoid,
                                                     unsigned partIndex) {
    std::vector<const double *> pcla;
    std::vector<const Node *> c = fromNode->GetChildren();
    for (auto i : c) {
        if (i != avoid) {
            InternalNodeWork * iw = (InternalNodeWork *)i->GetWork(partIndex);
            double * ic = &(iw->claAtParFromNd[0]);
            pcla.push_back(const_cast<const double *>(ic));
        }
    }
    if (fromNode->parent && fromNode->parent != avoid) {
        InternalNodeWork * iw = (InternalNodeWork *)fromNode->GetWork(partIndex);
        double * ic = &(iw->claAtNdFromPar[0]);
        pcla.push_back(const_cast<const double *>(ic));
    }
    return pcla;
}

inline std::vector<const double *> Arc::GetPrevCLAs(unsigned partIndex) const {
    std::vector<const double *> pcla;
    if (fromIsChild) {
        assert(fromNode->leftChild);
        std::vector<const Node *> c = const_cast<const Node *>(fromNode)->GetChildren();
        for (auto i : c) {
            InternalNodeWork * iw = (InternalNodeWork *)i->GetWork(partIndex);
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
