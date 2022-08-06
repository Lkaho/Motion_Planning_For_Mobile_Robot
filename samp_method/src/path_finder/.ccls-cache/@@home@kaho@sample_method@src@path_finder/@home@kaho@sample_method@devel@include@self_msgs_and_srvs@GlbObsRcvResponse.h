// Generated by gencpp from file self_msgs_and_srvs/GlbObsRcvResponse.msg
// DO NOT EDIT!


#ifndef SELF_MSGS_AND_SRVS_MESSAGE_GLBOBSRCVRESPONSE_H
#define SELF_MSGS_AND_SRVS_MESSAGE_GLBOBSRCVRESPONSE_H


#include <string>
#include <vector>
#include <map>

#include <ros/types.h>
#include <ros/serialization.h>
#include <ros/builtin_message_traits.h>
#include <ros/message_operations.h>


namespace self_msgs_and_srvs
{
template <class ContainerAllocator>
struct GlbObsRcvResponse_
{
  typedef GlbObsRcvResponse_<ContainerAllocator> Type;

  GlbObsRcvResponse_()
    : result(false)  {
    }
  GlbObsRcvResponse_(const ContainerAllocator& _alloc)
    : result(false)  {
  (void)_alloc;
    }



   typedef uint8_t _result_type;
  _result_type result;





  typedef boost::shared_ptr< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> > Ptr;
  typedef boost::shared_ptr< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> const> ConstPtr;

}; // struct GlbObsRcvResponse_

typedef ::self_msgs_and_srvs::GlbObsRcvResponse_<std::allocator<void> > GlbObsRcvResponse;

typedef boost::shared_ptr< ::self_msgs_and_srvs::GlbObsRcvResponse > GlbObsRcvResponsePtr;
typedef boost::shared_ptr< ::self_msgs_and_srvs::GlbObsRcvResponse const> GlbObsRcvResponseConstPtr;

// constants requiring out of line definition



template<typename ContainerAllocator>
std::ostream& operator<<(std::ostream& s, const ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> & v)
{
ros::message_operations::Printer< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> >::stream(s, "", v);
return s;
}


template<typename ContainerAllocator1, typename ContainerAllocator2>
bool operator==(const ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator1> & lhs, const ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator2> & rhs)
{
  return lhs.result == rhs.result;
}

template<typename ContainerAllocator1, typename ContainerAllocator2>
bool operator!=(const ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator1> & lhs, const ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator2> & rhs)
{
  return !(lhs == rhs);
}


} // namespace self_msgs_and_srvs

namespace ros
{
namespace message_traits
{





template <class ContainerAllocator>
struct IsMessage< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> >
  : TrueType
  { };

template <class ContainerAllocator>
struct IsMessage< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> const>
  : TrueType
  { };

template <class ContainerAllocator>
struct IsFixedSize< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> >
  : TrueType
  { };

template <class ContainerAllocator>
struct IsFixedSize< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> const>
  : TrueType
  { };

template <class ContainerAllocator>
struct HasHeader< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> >
  : FalseType
  { };

template <class ContainerAllocator>
struct HasHeader< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> const>
  : FalseType
  { };


template<class ContainerAllocator>
struct MD5Sum< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> >
{
  static const char* value()
  {
    return "eb13ac1f1354ccecb7941ee8fa2192e8";
  }

  static const char* value(const ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator>&) { return value(); }
  static const uint64_t static_value1 = 0xeb13ac1f1354ccecULL;
  static const uint64_t static_value2 = 0xb7941ee8fa2192e8ULL;
};

template<class ContainerAllocator>
struct DataType< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> >
{
  static const char* value()
  {
    return "self_msgs_and_srvs/GlbObsRcvResponse";
  }

  static const char* value(const ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator>&) { return value(); }
};

template<class ContainerAllocator>
struct Definition< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> >
{
  static const char* value()
  {
    return "bool result\n"
"\n"
;
  }

  static const char* value(const ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator>&) { return value(); }
};

} // namespace message_traits
} // namespace ros

namespace ros
{
namespace serialization
{

  template<class ContainerAllocator> struct Serializer< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> >
  {
    template<typename Stream, typename T> inline static void allInOne(Stream& stream, T m)
    {
      stream.next(m.result);
    }

    ROS_DECLARE_ALLINONE_SERIALIZER
  }; // struct GlbObsRcvResponse_

} // namespace serialization
} // namespace ros

namespace ros
{
namespace message_operations
{

template<class ContainerAllocator>
struct Printer< ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator> >
{
  template<typename Stream> static void stream(Stream& s, const std::string& indent, const ::self_msgs_and_srvs::GlbObsRcvResponse_<ContainerAllocator>& v)
  {
    s << indent << "result: ";
    Printer<uint8_t>::stream(s, indent + "  ", v.result);
  }
};

} // namespace message_operations
} // namespace ros

#endif // SELF_MSGS_AND_SRVS_MESSAGE_GLBOBSRCVRESPONSE_H
