// file      : xsd/usage.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef USAGE_HXX
#define USAGE_HXX

#include <cult/types.hxx>

#include <backend-elements/indentation/buffer.hxx>
#include <backend-elements/indentation/clip.hxx>

namespace CLI
{
  using namespace Cult::Types;

  template <typename C>
  class OptionsUsage: public BackendElements::Indentation::Buffer<C>
  {
    typedef BackendElements::Indentation::Buffer<C> Buffer;

  public:
    typedef
    typename Buffer::Traits
    Traits;

    typedef
    typename Buffer::AsChar
    AsChar;

    typedef
    typename Buffer::AsInt
    AsInt;

    typedef
    typename Buffer::Write
    Write;

  public:
    OptionsUsage (Buffer& out)
        : out_ (out),
          option_length_ (0),
          construct_ (Construct::newline)
    {
    }

  public:
    virtual AsInt
    put (AsChar c)
    {
      AsInt result (Traits::to_int_type (c));

      try
      {
        switch (c)
        {
        case '\n':
          {
            switch (construct_)
            {
            case Construct::newline:
              {
                result = out_.put (c);
                break;
              }
            case Construct::option:
              {
                construct_ = Construct::newline;
                break;
              }
            case Construct::description:
              {
                result = out_.put (c);
                construct_ = Construct::newline;
                break;
              }
            default:
              {
                abort ();
              }
            }

            break;
          }
        case '-':
          {
            switch (construct_)
            {
            case Construct::newline:
              {
                construct_ = Construct::option;

                option_length_ = 0;

                output_indentation ();
                result = out_.put (c);

                ++option_length_;

                break;
              }
            case Construct::option:
              {
                ++option_length_;
                //fall through
              }
            case Construct::description:
              {
                result = out_.put (c);
                break;
              }
            default:
              {
                abort ();
              }
            }

            break;
          }
        default:
          {
            switch (construct_)
            {
            case Construct::newline:
              {
                construct_ = Construct::description;

                output_indentation ();

                result = out_.put (c);
                break;
              }
            case Construct::option:
              {
                ++option_length_;
                //fall through
              }
            default:
              {
                result = out_.put (c);
                break;
              }
            }

            break;
          }
        }
      }
      catch (Write const&)
      {
        result = Traits::eof ();
      }

      return result;
    }

    virtual Void
    unbuffer ()
    {
    }

  private:
    Void
    output_indentation ()
    {
      UnsignedLong spaces;

      switch (construct_)
      {
      case Construct::option:
        {
          spaces = 2;
          option_length_ += 2;
          break;
        }
      case Construct::description:
        {
          spaces = 29;

          if (option_length_)
          {
            if (option_length_ > spaces)
              spaces = 1;
            else
              spaces -= option_length_;

            option_length_ = 0;
          }

          break;
        }
      default:
        {
          abort ();
        }
      }

      while (spaces--)
        out_.put (' ');
    }

  private:
    Buffer& out_;
    UnsignedLong option_length_;

    struct Construct
    {
      enum Value
      {
        newline,
        option,
        description
      };
    };

    typename Construct::Value construct_;
  };

  //@@ rename Indentation to Indent in be?
  //
  namespace Indent = BackendElements::Indentation;
}

#endif // USAGE_HXX

