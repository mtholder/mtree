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
        double * GetCLA(unsigned i);

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
class NodeIterator {
    public:
        NodeIterator(Node *c) 
            : curr(c) {
        }
        virtual ~NodeIterator(){}
        Node * get() {
            return this->curr;
        }
        Node * next() {
            this->advance();
            return this->get();
        }
        virtual void advance() = 0;
    protected:
        Node * curr;
};
class PostorderNodeIterator:public NodeIterator {
    public:
        PostorderNodeIterator(Node * r, Node * avoidNode)
            :NodeIterator(r),
             avoid(avoidNode) {
            this->reset(r, avoidNode);
        }
        void reset(Node *r, Node *avoidNode) {
            while (!ancStack.empty()) {
                ancStack.pop();
            }
            avoid = avoidNode;
            curr = r;
            if (curr && curr->leftChild != nullptr) {
                if (curr->leftChild == avoid) {
                    ancStack.push(curr);
                    curr = curr->leftChild->rightSib;
                    if (curr == nullptr) {
                        curr = r;
                    } else {
                        this->addAncs(curr);
                    }
                } else {
                    this->addAncs(curr);
                }
            }
        }
        void advance() {
            if (curr->rightSib) {
                if (curr->rightSib == avoid) {
                    if (avoid->rightSib) {
                        curr = avoid->rightSib;
                        this->addAncs(curr);
                        return;
                    }
                } else {
                    curr = curr->rightSib;
                    this->addAncs(curr);
                    return;
                }
            }
            if (ancStack.empty()) {
                curr = nullptr;
            } else {
                curr = ancStack.top();
                ancStack.pop();
            }
        }
    private:
        void addAncs(Node *c) {
            curr = c;
            while (curr->leftChild) {
                ancStack.push(curr);
                curr = curr->leftChild;
            }
        }
    Node * avoid;
    std::stack<Node *>ancStack;

};
class PostorderForNodeIterator: public NodeIterator {
    public:
        PostorderForNodeIterator(Node * vr)
            :NodeIterator(vr),
            refNode(vr),
            post(nullptr, nullptr) {
            assert(vr);
            assert(vr->parent);
            curr = vr;
            while (curr->parent) {
                toRoot.push(curr);
                curr = curr->parent;
            }
            avoid = toRoot.top();
            toRoot.pop();
            post.reset(curr, avoid);
            curr = post.get();
            belowNode = true;
        }
        void advance() {
            curr = post.get();
            if (curr == nullptr && belowNode) {
                if (avoid == refNode) {
                    belowNode = false;
                    post.reset(refNode, nullptr);
                    curr = post.get();
                } else {
                    curr = avoid;
                    avoid = toRoot.top();
                    toRoot.pop();
                    post.reset(curr, avoid);
                    curr = post.get();
                }
                assert(curr != nullptr);
            }
        }
    private:
        Node * refNode;
        Node * avoid;
        std::stack<Node *> toRoot;
        PostorderNodeIterator post;
        bool belowNode;
};

inline NodeIterator * postorder(Node *c) {
    if (c->parent) {
        return new PostorderForNodeIterator(c);
    }
    return new PostorderNodeIterator(c, nullptr);
}

class LeafWork {
    public:
        LeafWork(unsigned numChars, unsigned numStateCodes, unsigned numStates, unsigned numRates) 
            :summed(numStateCodes*numStates*numRates),
             cla(numStates*numRates*numChars) {
        }
    std::vector<double> summed;
    std::vector<double> cla;
};

class InternalNodeWork {
    public:
        InternalNodeWork(unsigned numChars, unsigned numStates, unsigned numRates) 
            :cla(numChars*numStates*numRates),
            nChars(numChars) {
        }
    std::vector<double> cla;
    unsigned nChars;
};
inline double * Node::GetCLA(unsigned i) {
    void * w = this->GetWork(i);
    if (this->IsLeaf()) {
        LeafWork * lw = (LeafWork *)(w);
        return &(lw->cla[0]);
    } else {
        InternalNodeWork * iw = (InternalNodeWork*)(w);
        return &(iw->cla[0]);
    }
}
}
#endif
