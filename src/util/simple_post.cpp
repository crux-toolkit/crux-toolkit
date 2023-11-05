//#define BOOST_ASIO_ENABLE_HANDLER_TRACKING
#define BOOST_BIND_GLOBAL_PLACEHOLDERS

#include <iostream>
#include <iomanip>
#include <boost/bind.hpp>
#include <boost/asio.hpp>
#include <boost/asio/ssl.hpp>

class client
{
public:
  client(
      boost::asio::io_service& io_service,
      boost::asio::ssl::context& context,
      boost::asio::ip::tcp::resolver::iterator endpoint_iterator
  ) : socket_(io_service, context)
  {
    boost::asio::async_connect(socket_.lowest_layer(), endpoint_iterator,
        boost::bind(&client::handle_connect, this,
          boost::asio::placeholders::error));
  }

  void handle_connect(const boost::system::error_code& error)
  {
      if (!error)
      {
          socket_.async_handshake(boost::asio::ssl::stream_base::client,
                  boost::bind(&client::handle_handshake, this,
                      boost::asio::placeholders::error));
      }
      else
      {
          std::cout << "Connect failed: " << error.message() << "\n";
      }
  }

  void handle_handshake(const boost::system::error_code& error)
  {
      if (!error)
      {
          // Construct JSON for POST
          std::stringstream body;
          body
            << "{"
            << "  \"client_id\" : \"332557735.1693348426\","
            << "  \"events\" : ["
            << "    {"
            << "      \"name\" : \"crux\","
            << "      \"params\" : {"
            << "        \"tool\" : \"comet\","
            << "        \"platform\" : \"linux\","
            << "        \"version\" : \"4.1-bc4562ef-2023-07-25\""
            << "      }"
            << "    }"
            << "  ]"
            << "}";
          body.seekp(0, std::ios::end);
          std::stringstream::pos_type body_length = body.tellp();

          // Contruct POST headers and append JSON as body
          std::stringstream post_content(""); 
          post_content
            << "POST /mp/collect?"
            << "measurement_id=G-V7XKGGFPYX"
            << "&api_secret=UIf4l54KSbK84hPRWng2Yg"
            << " HTTP/1.1\n"
            << "Host: www.google-analytics.com\n"
            << "accept: */*\n"
            << "Content-Type: application/json\n"
            << "Content-Length: " << body_length
            << "\n\n"
            << body.str();

          // Request the POST
          boost::asio::async_write(socket_,
                  boost::asio::buffer(post_content.str()),
                  boost::bind(&client::handle_write, this,
                      boost::asio::placeholders::error,
                      boost::asio::placeholders::bytes_transferred));

      }
      else
      {
          std::cout << "Handshake failed: " << error.message() << "\n";
      }
  }

  void handle_write(const boost::system::error_code& error,
      size_t /*bytes_transferred*/)
  {
      // std::cerr << "starting read loop\n";
      if (!error)
      {
          boost::asio::async_read_until(socket_,
                  //boost::asio::buffer(reply_, sizeof(reply_)),
                  reply_, '\n',
                  boost::bind(&client::handle_read, this,
                      boost::asio::placeholders::error,
                      boost::asio::placeholders::bytes_transferred));
      }
      else
      {
          std::cout << "Write failed: " << error.message() << "\n";
      }
  }

  void handle_read(const boost::system::error_code& error, size_t /*bytes_transferred*/)
  {
      if (!error)
      {
          std::istream is(&reply_);
          std::string line;
          std::getline(is, line);
          //std::cout << line << "\n";
      }
      else
      {
          std::cout << "Read failed: " << error.message() << "\n";
      }
  }

private:
  boost::asio::ssl::stream<boost::asio::ip::tcp::socket> socket_;
  boost::asio::streambuf reply_;
};

int main(int argc, char* argv[])
{
    try
    {
        boost::asio::io_service io_service;

        boost::asio::ip::tcp::resolver resolver(io_service);
        std::string host = "www.google-analytics.com";
        std::string protocol = "https";
        boost::asio::ip::tcp::resolver::query query(host.c_str(), protocol.c_str());
        boost::asio::ip::tcp::resolver::iterator iterator = resolver.resolve(query);
        boost::asio::ssl::context ctx(boost::asio::ssl::context::sslv23);
        ctx.set_default_verify_paths();

        client c(io_service, ctx, iterator);

        io_service.run();
    }
    catch (std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << "\n";
    }

    return 0;
}
