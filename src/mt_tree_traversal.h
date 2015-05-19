#if !defined(__TREE_TRAVERSAL_H__)
#define __TREE_TRAVERSAL_H__
#include "mt_tree.h"
namespace mt {

class ArcIterator : std::forward_iterator_tag {
    public:
        ArcIterator(Node *c) 
            : currNd(c) {
        }
        virtual ~ArcIterator(){}
        virtual Arc get() const = 0;
        Arc next() {
            this->_advance();
            return this->get();
        }
        virtual void _advance() = 0;
        Node * GetCurrNode() {
            return currNd;
        }
        Node * currNd;

};

class PostorderArcIterator:public ArcIterator {
    public:
        PostorderArcIterator(Node * r, Node * avoidNode=nullptr)
            :ArcIterator(r),
             avoid(avoidNode) {
            if (r != nullptr) {
                this->reset(r, avoidNode);
            } else {
                this->currNd = avoidNode;
            }
        }
        bool operator==(const PostorderArcIterator &other) const {
            return ((this->currNd == other.currNd)
                    && (this->ancStack == other.ancStack)
                    && (this->avoid == other.avoid));
        }
        bool operator!=(const PostorderArcIterator &other) const {
            return ! (*this == other);
        }
        Arc operator*() const {
            return get();
        }
        PostorderArcIterator & operator++() {
            if (currNd == nullptr && ancStack.empty()) {
                throw std::out_of_range("Incremented a dead PostorderArcIterator");
            }
            _advance();
            return *this;
        }
        virtual Arc get() const {
            if (ancStack.empty()) {
                return Arc(nullptr, this->currNd);
            }
            return Arc(this->currNd, this->currNd->parent);
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
        void _advance() {
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
        Arc get() const {
            if (currNd == rootOfCurrPost) {
                return Arc(currNd, (currNd == refNode ? nullptr : avoid));
            } else {
                return post.get();
            }
        }
        void _advance() {
            post._advance();
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

class ConstPostorderForNodeIterator {
    public:
        ConstPostorderForNodeIterator(const Node * vr)
            :pfni(const_cast<Node*>(vr)){
        }
        ConstArc get() {
            return ConstArc(pfni.get());
        }
        ConstArc next() {
            return ConstArc{pfni.next()};
        }
        const Node * GetCurrNode() {
            return const_cast<const Node*>(pfni.GetCurrNode());
        }
    private:
        PostorderForNodeIterator pfni;

};

class PostArcIter {
    private:
    Node * startNode;
    public:
    explicit PostArcIter(Node *s)
        :startNode(s) {
    }
    PostorderArcIterator begin() const {
        return PostorderArcIterator{startNode, startNode};
    }
    PostorderArcIterator end() const {
        return PostorderArcIterator{nullptr, startNode};
    }
};


inline PostorderForNodeIterator postorder(Node *c) {
    return PostorderForNodeIterator(c);
}
inline ConstPostorderForNodeIterator postorder(const Node *c) {
    return ConstPostorderForNodeIterator(c);
}

} // namespace mt
#endif
