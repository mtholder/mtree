#if !defined(__MT_INSTANCE_H__)
#   error "This impl file is included by mt_instance, and should not be included eleswhere."
#endif
#if !defined(__MT_INSTANCE_INLINES_H__)
#define __MT_INSTANCE_INLINES_H__
// This poorly-named file has inlines that requre MTInstance be included.
namespace mt {

inline SlotIndices::SlotIndices(const MTInstance &instance, const TraversalInfo & tInfo) 
    :p(instance.GetUseRecom() ? tInfo.slot_p : tInfo.pNumber - instance.GetMxTips() - 1),
    q(instance.GetUseRecom() ? tInfo.slot_p : tInfo.pNumber - instance.GetMxTips() - 1),
    r(instance.GetUseRecom() ? tInfo.slot_p : tInfo.pNumber - instance.GetMxTips() - 1) {
}

inline bool PartModelInfo::GetUseRecom() const{
    return instance.GetUseRecom();
}

} // namespace
#endif
