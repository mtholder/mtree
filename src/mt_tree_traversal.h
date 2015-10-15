#if !defined(__TREE_TRAVERSAL_H__)
#define __TREE_TRAVERSAL_H__
#include "mt_tree.h"
namespace mt {

double ScoreTree(PartitionedMatrix &partMat, Tree &tree, CharModel &cm);

class ArcIterator {
    public:
        ArcIterator()
            : currNd(nullptr) {
        }
        ArcIterator(Node *c)
            : currNd(c) {
        }
        virtual ~ArcIterator(){}
        virtual Arc get() = 0;
        Arc next() {
            this->advance();
            return this->get();
        }
        virtual void advance() = 0;
        Node * GetCurrNode() {
            return currNd;
        }
        Node * currNd = nullptr;
};

class PostorderArcIterator:public ArcIterator {
    public:
        PostorderArcIterator(Node * r, Node * avoidNode)
            :ArcIterator(r),
             avoid(avoidNode) {
            this->reset(r, avoidNode);
        }
        virtual Arc get() {
            return Arc(this->currNd,
                       (ancStack.empty() ? nullptr : this->currNd->parent));
        }
        void reset(Node *r, Node *avoidNode) {
            while (!ancStack.empty()) {
                ancStack.pop();
            }
            avoid = avoidNode;
            currNd = r;
            if (currNd && currNd->leftChild != nullptr) {
                if (currNd->leftChild == avoid) {
                    ancStack.push(currNd);
                    currNd = currNd->leftChild->rightSib;
                    if (currNd == nullptr) {
                        currNd = r;
                    } else {
                        this->addAncs(currNd);
                    }
                } else {
                    this->addAncs(currNd);
                }
            }
        }
        void advance() {
            if (ancStack.empty()) {
                currNd = nullptr;
                return;
            }
            if (currNd->rightSib) {
                if (currNd->rightSib == avoid) {
                    if (avoid->rightSib) {
                        currNd = avoid->rightSib;
                        this->addAncs(currNd);
                        return;
                    }
                } else {
                    currNd = currNd->rightSib;
                    this->addAncs(currNd);
                    return;
                }
            }
            currNd = ancStack.top();
            ancStack.pop();
        }
    private:
        void addAncs(Node *c) {
            currNd = c;
            while (currNd->leftChild) {
                ancStack.push(currNd);
                currNd = currNd->leftChild;
            }
        }
        Node * avoid;
        std::stack<Node *>ancStack;

};


class PostorderForNodeIterator: public ArcIterator {
    public:
        PostorderForNodeIterator(Node * vr)
            :ArcIterator(vr),
            refNode(vr),
            rootOfCurrPost(nullptr),
            post(nullptr, nullptr) {
            assert(vr);
            currNd = vr;
            if (currNd->parent) {
                while (currNd->parent) {
                    toRoot.push(currNd);
                    currNd = currNd->parent;
                }
                avoid = toRoot.top();
                toRoot.pop();
                post.reset(currNd, avoid);
                rootOfCurrPost = currNd;
                currNd = post.GetCurrNode();
                belowNode = true;
            } else {
                belowNode = false;
                post.reset(refNode, nullptr);
                rootOfCurrPost = refNode;
                currNd = post.GetCurrNode();
            }
        }
        Arc get() {
            if (currNd == rootOfCurrPost) {
                return Arc(currNd, (currNd == refNode ? nullptr : avoid));
            } else {
                return post.get();
            }
        }
        void advance() {
            post.advance();
            currNd = post.GetCurrNode();
            if (currNd == nullptr && belowNode) {
                if (avoid == refNode) {
                    belowNode = false;
                    post.reset(refNode, nullptr);
                    rootOfCurrPost = refNode;
                    currNd = post.GetCurrNode();
                } else {
                    currNd = avoid;
                    avoid = toRoot.top();
                    toRoot.pop();
                    post.reset(currNd, avoid);
                    rootOfCurrPost = currNd;
                    currNd = post.GetCurrNode();
                }
                assert(currNd != nullptr);
            }
        }
    private:
        Node * refNode;
        Node * avoid;
        Node * rootOfCurrPost;
        std::stack<Node *> toRoot;
        PostorderArcIterator post;
        bool belowNode;
};

inline PostorderForNodeIterator postorder(Node *c) {
    return PostorderForNodeIterator(c);
}
}
#endif
