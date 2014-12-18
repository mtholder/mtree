#if !defined(__TREE_H__)
#define __TREE_H__
#include <cassert>
#include <stack>
#include <climits>
#include <iostream>
#include <vector>
namespace mt {

typedef unsigned int char_state_t;
class CharStateToPrimitiveInd {
    public:
        CharStateToPrimitiveInd(unsigned numStateCodes)
            :stateCodeToStateCodeVec(numStateCodes) {
        }
        const std::vector<char_state_t> & GetStateCodes(unsigned i) const {
            return stateCodeToStateCodeVec[i];
        }
        void SetStateCode(unsigned i, const std::vector<char_state_t> & v) {
            this->stateCodeToStateCodeVec[i] = v;
        }
    private:
        std::vector<std::vector<char_state_t> > stateCodeToStateCodeVec;
};
class LeafCharacterVector {
    public:
        LeafCharacterVector()
            :cs2pi(nullptr) {
        }
        LeafCharacterVector(const char_state_t *inp, unsigned len, const CharStateToPrimitiveInd * stateToPrimStates)
            :charVec(len),
            cs2pi(stateToPrimStates) {
            for (auto i = 0U; i < len; ++i) {
                charVec[i] = inp[i];
            }
        }
        std::vector<char_state_t> charVec;
        const CharStateToPrimitiveInd * cs2pi;
};

typedef std::vector<LeafCharacterVector> CharMatrix;

class PartitionedMatrix {
    public:
        PartitionedMatrix(unsigned numTaxa, const std::vector<unsigned> &numCharsPerPartition, const std::vector<unsigned> &orig2compressed)
            :nTaxa(numTaxa),
             nCharsVec(numCharsPerPartition),
             partitions(numCharsPerPartition.size()),
             origInd2CompressedInd(orig2compressed) {
             for (auto i : partitions) {
                i.resize(numTaxa);
             }
        }
        unsigned GetNumPartitions() const {
            return this->partitions.size();
        }
        const LeafCharacterVector * GetLeafCharacters(unsigned partIndex, unsigned leafIndex) const {
            return &(this->partitions[partIndex][leafIndex]);
        }
        void fillPartition(unsigned partIndex, const char_state_t ** mat, const CharStateToPrimitiveInd *cs2pi) {
            CharMatrix & part = partitions[partIndex];
            part.resize(nTaxa);
            for (auto i = 0U; i < nTaxa; ++i) {
                part[i] = LeafCharacterVector(mat[i], nCharsVec[partIndex], cs2pi);
            }
        }
    private:
        unsigned nTaxa;
        std::vector<unsigned> nCharsVec;
        std::vector<CharMatrix> partitions;
        std::vector<unsigned> origInd2CompressedInd;

};
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
        void * GetWork(unsigned i) {
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
        unsigned GetsNumLeaves() const {
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

class CharModel {
    public:
        CharModel(unsigned numStates, unsigned numRateCats)
            :nStates(numStates),
            nRateCats(numRateCats) {
        }
        virtual ~CharModel() {
        }
        virtual double sumLnL(const Node * virtualRoot) const = 0;
        virtual void fillLeafWork(const LeafCharacterVector *, LeafWork *, double edgeLen);
        virtual void fillInternalWork(const double * cla1, const double *cla2, InternalNodeWork *, double edgeLen);
        virtual double * calcTransitionProb(double edgeLen) = 0;
    protected:
        unsigned nStates;
        unsigned nRateCats;
};
class MkVarNoMissingAscCharModel: public CharModel {
    public:
        MkVarNoMissingAscCharModel(unsigned numStates, unsigned numRateCats)
            :CharModel(numStates, numRateCats) {
        }
        virtual ~MkVarNoMissingAscCharModel() {
        }
        virtual double sumLnL(const Node * ) const {
            return 3.2;
        }
        virtual double * calcTransitionProb(double edgeLen) {
            return nullptr;
        }
};
class MkCharModel: public CharModel {
    public:
        MkCharModel(unsigned numStates, unsigned numRateCats)
            :CharModel(numStates, numRateCats) {
        }
        virtual ~MkCharModel() {
        }
        virtual double sumLnL(const Node * ) const {
            return 3.4;
        }
        virtual double * calcTransitionProb(double edgeLen) {
            return nullptr;
        }
};
void doAnalysis(Tree &tree, CharModel &cm);

class InternalNodeWork {
    public:
        InternalNodeWork(unsigned numChars, unsigned numStates, unsigned numRates) 
            :claAtNdFromNd(numStates*numRates*numChars),
             claAtNdFromPar(numStates*numRates*numChars),
             claAtParFromNd(numStates*numRates*numChars),
             claAtParFromPar(numStates*numRates*numChars) {
        }
        unsigned GetLenCLA() const {
            return claAtNdFromNd.size();
        }
    std::vector<double> claAtNdFromNd;
    std::vector<double> claAtNdFromPar;
    std::vector<double> claAtParFromNd;
    std::vector<double> claAtParFromPar;
    unsigned nChars;
};

class LeafWork: public InternalNodeWork {
    public:
        LeafWork(unsigned numChars, unsigned numStateCodes, unsigned numStates, unsigned numRates) 
            :InternalNodeWork(numChars, numStates, numRates),
            summed(numStateCodes*numStates*numRates){}
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
            return (const LeafCharacterVector *)fromNode->GetData(partIndex);
        }
        LeafWork * GetFromNdLeafWork(unsigned partIndex) {
            assert(IsFromLeaf());
            return (LeafWork *)fromNode->GetWork(partIndex);
        }
        InternalNodeWork * GetFromNdIntWork(unsigned partIndex) {
            return (InternalNodeWork *)toNode->GetWork(partIndex);
        }
        InternalNodeWork * GetToNdIntWork(unsigned partIndex) {
            return (InternalNodeWork *)toNode->GetWork(partIndex);
        }
        double * GetFromNdCLA(unsigned partIndex, bool crossedEdge);
        double * GetToNdCLA(unsigned partIndex, bool crossedEdge);
        std::vector<const double *> GetPrevCLAs(unsigned partIndex);
        unsigned GetLenCLA(unsigned partIndex) {
            return GetToNdIntWork(partIndex)->GetLenCLA();
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

inline std::vector<const double *> Arc::GetPrevCLAs(unsigned partIndex) {
    std::vector<const double *> pcla;
    if (fromIsChild) {
        assert(fromNode->leftChild);
        std::vector<Node *> c = fromNode->GetChildren();
        for (auto i : c) {
            InternalNodeWork * iw = (InternalNodeWork *)i->GetWork(partIndex);
            double * ic = &(iw->claAtParFromNd[0]);
            pcla.push_back(const_cast<const double *>(ic));
        }
    } else {
        std::vector<Node *> c = fromNode->GetChildren();
        for (auto i : c) {
            if (i != toNode) {
                InternalNodeWork * iw = (InternalNodeWork *)i->GetWork(partIndex);
                double * ic = &(iw->claAtParFromNd[0]);
                pcla.push_back(const_cast<const double *>(ic));
            }
        }
        if (fromNode->parent) {
            InternalNodeWork * iw = (InternalNodeWork *)fromNode->GetWork(partIndex);
            double * ic = &(iw->claAtNdFromPar[0]);
            pcla.push_back(const_cast<const double *>(ic));
        }
    }
    return pcla;
}


} //namespace
#endif
