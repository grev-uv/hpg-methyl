#ifndef _CP_HTTP_H
#define _CP_HTTP_H

/**
 * @addtogroup cp_socket
 */
/** @{ */
/**
 * @file
 * basic http feature implementations. the intention is to offer a simple way
 * of providing services to http based clients. 
 */

#include "config.h"
#include "common.h"

#ifdef _WINDOWS
#include <Winsock2.h>
#endif

__BEGIN_DECLS

#include "hashtable.h"
#include "trie.h"
#include "vector.h"
#include "socket.h"
#include "str.h"

#include <string.h>
#include <time.h>

/* HTTP 1.1 defines these methods (rfc 2616)
 * 
 *     Method         = "OPTIONS"
 *                    | "GET"   
 *                    | "HEAD"  
 *                    | "POST"  
 *                    | "PUT"   
 *                    | "DELETE"
 *                    | "TRACE" 
 *                    | "CONNECT"
 *
 *     general-header = Cache-Control     
 *                    | Connection        
 *                    | Date              
 *                    | Pragma            
 *                    | Trailer           
 *                    | Transfer-Encoding 
 *                    | Upgrade           
 *                    | Via               
 *                    | Warning           
 *
 *
 * 
 *     request-header = Accept             
 *                    | Accept-Charset     
 *                    | Accept-Encoding    
 *                    | Accept-Language    
 *                    | Authorization      
 *                    | Expect             
 *                    | From               
 *                    | Host               
 *                    | If-Match           
 *                    | If-Modified-Since  
 *                    | If-None-Match      
 *                    | If-Range           
 *                    | If-Unmodified-Since
 *                    | Max-Forwards       
 *                    | Proxy-Authorization
 *                    | Range              
 *                    | Referer            
 *                    | TE                 
 *                    | User-Agent         
 *
 *
 *     entity-header  = Allow              
 *                    | Content-Encoding   
 *                    | Content-Language   
 *                    | Content-Length     
 *                    | Content-Location   
 *                    | Content-MD5        
 *                    | Content-Range      
 *                    | Content-Type       
 *                    | Expires            
 *                    | Last-Modified      
 *
 *
 *     Status-Code    = "100"  : Continue
 *                    | "101"  : Switching Protocols
 *                    | "200"  : OK
 *                    | "201"  : Created
 *                    | "202"  : Accepted
 *                    | "203"  : Non-Authoritative Information
 *                    | "204"  : No Content
 *                    | "205"  : Reset Content
 *                    | "206"  : Partial Content
 *                    | "300"  : Multiple Choices
 *                    | "301"  : Moved Permanently
 *                    | "302"  : Found
 *                    | "303"  : See Other
 *                    | "304"  : Not Modified
 *                    | "305"  : Use Proxy
 *                    | "307"  : Temporary Redirect
 *                    | "400"  : Bad Request
 *                    | "401"  : Unauthorized
 *                    | "402"  : Payment Required
 *                    | "403"  : Forbidden
 *                    | "404"  : Not Found
 *                    | "405"  : Method Not Allowed
 *                    | "406"  : Not Acceptable
 *                    | "407"  : Proxy Authentication Required
 *                    | "408"  : Request Time-out
 *                    | "409"  : Conflict
 *                    | "410"  : Gone
 *                    | "411"  : Length Required
 *                    | "412"  : Precondition Failed
 *                    | "413"  : Request Entity Too Large
 *                    | "414"  : Request-URI Too Large
 *                    | "415"  : Unsupported Media Type
 *                    | "416"  : Requested range not satisfiable
 *                    | "417"  : Expectation Failed
 *                    | "500"  : Internal Server Error
 *                    | "501"  : Not Implemented
 *                    | "502"  : Bad Gateway
 *                    | "503"  : Service Unavailable
 *                    | "504"  : Gateway Time-out
 *                    | "505"  : HTTP Version not supported
 *
 * 
 *
 *    response-header = Accept-Ranges           ; Section 14.5
 *                    | Age                     ; Section 14.6
 *                    | ETag                    ; Section 14.19
 *                    | Location                ; Section 14.30
 *                    | Proxy-Authenticate      ; Section 14.33
 *                    | Retry-After             ; Section 14.37
 *                    | Server                  ; Section 14.38
 *                    | Vary                    ; Section 14.44
 *                    | WWW-Authenticate        ; Section 14.47
 *
 */

/** HTTP request types */
typedef CPROPS_DLL enum 
{ 
	OPTIONS, 
	GET, 
	HEAD, 
	POST, 
	PUT, 
#ifndef _WINDOWS
	DELETE, 
#endif
	TRACE, 
	CONNECT 
} cp_http_request_type;

CPROPS_DLL
char *get_http_request_type_lit(cp_http_request_type type);

/** HTTP versions */
typedef CPROPS_DLL enum { HTTP_1_0, HTTP_1_1} cp_http_version;

/** HTTP status codes */
typedef CPROPS_DLL enum 
{
	HTTP_NULL_STATUS = -1,
	HTTP_100_CONTINUE = 100,
	HTTP_101_SWITCHING_PROTOCOLS = 101,
	HTTP_200_OK = 200,
	HTTP_201_CREATED = 201,
	HTTP_202_ACCEPTED = 202,
	HTTP_203_NON_AUTHORITATIVE_INFORMATION = 203,
	HTTP_204_NO_CONTENT = 204,
	HTTP_205_RESET_CONTENT = 205,
	HTTP_206_PARTIAL_CONTENT = 206,
	HTTP_300_MULTIPLE_CHOICES = 300,
	HTTP_301_MOVED_PERMANENTLY = 301,
	HTTP_302_FOUND = 302,
	HTTP_303_SEE_OTHER = 303,
	HTTP_304_NOT_MODIFIED = 304,
	HTTP_305_USE_PROXY = 305,
	HTTP_307_TEMPORARY_REDIRECT = 307,
	HTTP_400_BAD_REQUEST = 400,
	HTTP_401_UNAUTHORIZED = 401,
	HTTP_402_PAYMENT_REQUIRED = 402,
	HTTP_403_FORBIDDEN = 403,
	HTTP_404_NOT_FOUND = 404,
	HTTP_405_METHOD_NOT_ALLOWED = 405,
	HTTP_406_NOT_ACCEPTABLE = 406,
	HTTP_407_PROXY_AUTHENTICATION_REQUIRED = 407,
	HTTP_408_REQUEST_TIME_OUT = 408,
	HTTP_409_CONFLICT = 409,
	HTTP_410_GONE = 410,
	HTTP_411_LENGTH_REQUIRED = 411,
	HTTP_412_PRECONDITION_FAILED = 412,
	HTTP_413_REQUEST_ENTITY_TOO_LARGE = 413,
	HTTP_414_REQUEST_URI_TOO_LARGE = 414,
	HTTP_415_UNSUPPORTED_MEDIA_TYPE = 415,
	HTTP_416_REQUESTED_RANGE_NOT_SATISFIABLE = 416,
	HTTP_417_EXPECTATION_FAILED = 417,
	HTTP_500_INTERNAL_SERVER_ERROR = 500,
	HTTP_501_NOT_IMPLEMENTED = 501,
	HTTP_502_BAD_GATEWAY = 502,
	HTTP_503_SERVICE_UNAVAILABLE = 503,
	HTTP_504_GATEWAY_TIME_OUT = 504,
	HTTP_505_HTTP_VERSION_NOT_SUPPORTED = 505
} cp_http_status_code;

/** 
 * connection handling policy. HTTP_CONNECTION_POLICY_DEFAULT defaults to 
 * HTTP_CONNECTION_POLICY_CLOSE for HTTP/1.0 and to 
 * HTTP_CONNECTION_POLICY_KEEP_ALIVE for HTTP/1.1.
 */
typedef CPROPS_DLL enum
{ 
	HTTP_CONNECTION_POLICY_DEFAULT,
	HTTP_CONNECTION_POLICY_CLOSE, 
  	HTTP_CONNECTION_POLICY_KEEP_ALIVE
} connection_policy;

typedef CPROPS_DLL enum
{
	TEXT,
	HTML,
	JPEG
} cp_http_content_type;

#define DEFAULT_SERVER_NAME "cprops-0.1.12"
#define DEFAULT_KEEPALIVE 300

#ifndef HTTP_KEEPALIVE
#define HTTP_KEEPALIVE DEFAULT_KEEPALIVE
#endif

#ifndef MAX_URL_LENGTH
#define MAX_URL_LENGTH   0x400
#endif

/** recommended to call before using cp_httpsocket api */
CPROPS_DLL
int cp_http_init();

/** call to perform cleanup after using cp_httpsocket api */
CPROPS_DLL
void cp_http_shutdown();

#ifdef CP_USE_COOKIES

#ifndef MAX_COOKIE_LENGTH
#define MAX_COOKIE_LENGTH   0x1000
#endif /* MAX_COOKIE_LENGTH */

/* Wdy, DD-Mon-YYYY HH:MM:SS GMT */
#define COOKIE_TIME_FMT "a, 0-b-Y H:M: GMT"
#endif /* CP_USE_COOKIES */

#ifdef CP_USE_HTTP_SESSIONS

#define CP_HTTP_SESSION_PRM "CPSID"
#define CP_HTTP_SESSION_MARKER "CPSID="
#define CP_HTTP_SESSION_MARKER_LEN 6

#define DEFAULT_SESSION_VALIDITY 86400 /* 24 hours */

/** 
 * call from signal handler to stop sockets in waiting select() and close all 
 * connections
 */
CPROPS_DLL
void cp_httpsocket_stop_all();

/** 
 * cp_http_sessions are implemented with cookies if configured to use them and
 * supported by client or with url rewriting otherwise. 
 */
typedef CPROPS_DLL enum 
{ 
#ifdef CP_USE_COOKIES
	SESSION_TYPE_COOKIE = 1, 
#endif /* CP_USE_COOKIES */
	SESSION_TYPE_URLREWRITE = 2 
} cp_http_session_type;

typedef CPROPS_DLL struct _cp_http_session
{
	char *sid;					/**< reasonably unique session id          */
	cp_http_session_type type;	/**< session handling scheme               */
	time_t created;				/**< creation time                         */
	time_t access;				/**< last access                           */
	long validity;				/**< validity period in seconds            */
	short renew_on_access;      /**< validity calculated from last access? */
	cp_hashtable *key;			/**< session data table                    */
	short valid;				/**< currently valid                       */
	short fresh; 				/**< newly created                         */
	int refcount;				/**< reference count                       */
} cp_http_session;

/** retrieve the value stored on the session under the given key */
CPROPS_DLL
void *cp_http_session_get(cp_http_session *session, char *key);
/** set a (key, value) pair */
CPROPS_DLL
void *cp_http_session_set(cp_http_session *session, 
                                 char *key, void *value);
/** set a (key, value) pair with a given destructor function for the value */
CPROPS_DLL
void *cp_http_session_set_dtr(cp_http_session *session, 
		                      char *key, void *value, 
							  cp_destructor_fn dtr);
/** set validity period in seconds */
CPROPS_DLL
void cp_http_session_set_validity(cp_http_session *session, long sec);
/** session has not yet been joined by client */
CPROPS_DLL
short cp_http_session_is_fresh(cp_http_session *session);
/** (internal) deallocate an http session descriptor */
CPROPS_DLL
void cp_http_session_delete(cp_http_session *session);

#endif /* CP_USE_HTTP_SESSIONS */


struct CPROPS_DLL _cp_httpsocket;
		
/** http request descriptor */
typedef CPROPS_DLL struct _cp_http_request
{
	cp_http_request_type type; 		/**< GET, POST etc                     */
	cp_http_version version;		/**< HTTP_1_0, HTTP_1_1                */
	char *uri;						/**< request uri                       */
	cp_hashtable *header;			/**< table of headers                  */
	char *query_string;             /**< raw query string                  */
	cp_hashtable *prm;				/**< table of GET/POST parameters      */
	cp_vector *prm_name;			/**< parameter vector                  */
#ifdef CP_USE_COOKIES
	cp_vector *cookie;				/**< cookies                           */
#endif
	char *content;					/**< request content                   */
	struct _cp_httpsocket *owner;	/**< owning socket                     */
	cp_connection_descriptor *connection; /**< connection descriptor       */
#ifdef CP_USE_HTTP_SESSIONS
	cp_http_session *session;		/**< http session descriptor           */
#endif
} cp_http_request;

/** parse an http request */
CPROPS_DLL
cp_http_request *cp_http_request_parse(struct _cp_httpsocket *owner, 
		                               char *request, 
                                       int *err);
/** deallocate an http request descriptor */
CPROPS_DLL
void cp_http_request_delete(cp_http_request *request);
/** return value of specified header or NULL if none */
CPROPS_DLL
char *cp_http_request_get_header(cp_http_request *request, char *name);
/** return a newly allocated array of header names */
CPROPS_DLL
char **cp_http_request_get_headers(cp_http_request *request);
/** return a single parameter with the specified name */
CPROPS_DLL
char *cp_http_request_get_parameter(cp_http_request *request, char *name);
/** return vector of GET or POST parameters with specified name */
CPROPS_DLL
cp_vector *cp_http_request_get_param_vector(cp_http_request *request, char *name);
/** return all GET or POST parameters for given request */
CPROPS_DLL
cp_vector *cp_http_request_get_params(cp_http_request *request);
#ifdef CP_USE_HTTP_SESSIONS
/** return existing or create cp_http_session for given request */
CPROPS_DLL
cp_http_session *cp_http_request_get_session(cp_http_request *request, 
		                                     int create);
#endif
/** dump a cp_http_request descriptor (uses cp_info) */
CPROPS_DLL
void cp_http_request_dump(cp_http_request *req);

/** http response holder */
typedef CPROPS_DLL struct _cp_http_response
{
	cp_http_version version;			/**< HTTP_1_0, HTTP_1_1            */
	cp_http_status_code status;			/**< http return code (rfc 2616)   */
	cp_http_request *request;			/**< originating request           */
	char *servername;					/**< Server header                 */
	connection_policy connection;		/**< Keep-Alive, Close             */
	cp_hashtable *header;				/**< headers                       */
#ifdef CP_USE_COOKIES
	cp_vector *cookie;					/**< cookies                       */
#endif
	cp_http_content_type content_type; 	/**< deprecated                    */
	char *content_type_lit;				/**< content type string           */
	char *body; 						/**< deprecated                    */
	cp_string *content;					/**< content                       */
	int len;							
	short skip;							/**< skip flag                     */
	char *status_lit;                   /**< remote server code literal    */
} cp_http_response;

/** (internal) create a new response descriptor */
CPROPS_DLL
cp_http_response *cp_http_response_create(cp_http_request *request);
/** (internal) deallocate a response descriptor */
CPROPS_DLL
void cp_http_response_delete(cp_http_response *res);
/** deallocate a response descriptor */
CPROPS_DLL
void cp_http_response_destroy(cp_http_response *res);
/** (internal) output an http response on given connection */
CPROPS_DLL
int cp_http_response_write(cp_connection_descriptor *cdesc, 
						   cp_http_response *res);

/** set the status code on an http response */
CPROPS_DLL
void cp_http_response_set_status(cp_http_response *response, 
		                                cp_http_status_code code);
/** get the status code of an http response */
CPROPS_DLL
cp_http_status_code cp_http_response_get_status(cp_http_response *response);
/** deprecated - use cp_http_response_set_content_type_string */
CPROPS_DLL
void cp_http_response_set_content_type(cp_http_response *response, 
                                       cp_http_content_type type);
/** set the Content-Type header on an http response */
CPROPS_DLL
void cp_http_response_set_content_type_string(cp_http_response *response,
                                              char *content_type_lit);
/** get the Content-Type header value for an http response */
CPROPS_DLL
char *cp_http_response_get_content_type(cp_http_response *response);
/** set an arbitrary header on an http response */
CPROPS_DLL
void cp_http_response_set_header(cp_http_response *response, 
		                                char *name, char *value);
/** retrieve header content for one header */
CPROPS_DLL
char *cp_http_response_get_header(cp_http_response *response, char *name);
/** returns a vector containing header names */
CPROPS_DLL
cp_vector *cp_http_response_get_header_names(cp_http_response *response);
/** deprecated - use cp_http_response_set_content instead */
CPROPS_DLL
void cp_http_response_set_body(cp_http_response *response, 
	                                  char *body);
/** set the (possibly binary) content on an http response */
CPROPS_DLL
void cp_http_response_set_content(cp_http_response *response, 
		                                 cp_string *content);
/** get the content of an http response */
CPROPS_DLL
cp_string *cp_http_response_get_content(cp_http_response *response);

/** set the Connection header on a cp_http_response descriptor */
CPROPS_DLL
void cp_http_response_set_connection_policy(cp_http_response *response, 
		                                    connection_policy policy);

#ifdef CP_USE_COOKIES
/** add a Set-Cookie header to a cp_http_response descriptor */
CPROPS_DLL
int cp_http_response_set_cookie(cp_http_response *response, char *name, 
		char *content, char *host, char *path, long validity, int secure);
#endif

/** 
 * sets the 'skip' flag on an http response, which prevents the response being
 * deserialized and written to the client by the http framework. this is 
 * needed if the service function writes directly to the underlying socket.
 */
CPROPS_DLL
void cp_http_response_skip(cp_http_response *response);

/** set up a response reporting an error */
CPROPS_DLL
void cp_http_response_report_error(cp_http_response *response, 
		                           cp_http_status_code code, 
								   char *message);

/** service function prototype */
typedef int (*cp_http_service_callback)(cp_http_request *request, 
									    cp_http_response *response);

/** http server socket */
typedef CPROPS_DLL struct _cp_httpsocket
{
	int id;					/**< system assigned id                        */
	cp_socket *sock;		/**< underlying socket                         */
	int keepalive; 			/**< keep-alive in seconds                     */
#ifdef CP_USE_HTTP_SESSIONS
	cp_hashlist *session;	/**< active session table                      */
	cp_hashlist *pending;	/**< pending response list                     */
#endif
	cp_trie *service; 		/**< registered services (excluding default)   */
	cp_http_service_callback default_service;	/**< default service       */
	char *server_name;		/**< server name                               */
} cp_httpsocket;


/** 
 * http service descriptor
 */
typedef CPROPS_DLL struct _cp_http_service
{
	char *name; 						/**< descriptive name */
	char *path;							/**< base uri         */
	cp_http_service_callback service;   /**< service function */
} cp_http_service;

CPROPS_DLL
void *cp_http_thread_fn(void *prm);

/** create a cp_http_socket with the specified default cp_http_service */
CPROPS_DLL
cp_httpsocket *
    cp_httpsocket_create(int port, cp_http_service_callback default_service);
#ifdef CP_USE_SSL
/** 
 * create an ssl enabled cp_httpsocket structure. certificate_file and key_file
 * are paths to certificate and key in pem format.
 */
CPROPS_DLL
cp_httpsocket *
    cp_httpsocket_create_ssl(int port, 
			                 cp_http_service_callback default_service,
							 char *certificate_file, 
							 char *key_file, 
							 int verification_mode);
#endif
/** deallocate a cp_httpsocket structure */
CPROPS_DLL
void cp_httpsocket_delete(cp_httpsocket *svc);

/** set value for Keep-Alive header */
CPROPS_DLL
void cp_httpsocket_set_keepalive(cp_httpsocket *socket, int sec);
/** set value for Server header */
CPROPS_DLL
void cp_httpsocket_set_server_name(cp_httpsocket *socket, char *name);
/** 
 * set maximal number of queued requests before accept() starts refusing
 * connections
 */
CPROPS_DLL
void cp_httpsocket_set_backlog(cp_httpsocket *socket, int backlog);
/** set time to block in accept */
CPROPS_DLL
void cp_httpsocket_set_delay(cp_httpsocket *socket, struct timeval delay);
/* set sec. to block in accept() */
CPROPS_DLL
void cp_httpsocket_set_delay_sec(cp_httpsocket *socket, long sec);
/** set micro sec. to block in accept() */
CPROPS_DLL
void cp_httpsocket_set_delay_usec(cp_httpsocket *socket, long usec);
/** set lower limit on thread pool size */
CPROPS_DLL
void cp_httpsocket_set_poolsize_min(cp_httpsocket *socket, int min);
/** set upper limit on thread pool size */
CPROPS_DLL
void cp_httpsocket_set_poolsize_max(cp_httpsocket *socket, int max);

/** called when closing a socket */
CPROPS_DLL
void *cp_httpsocket_add_shutdown_callback(cp_httpsocket *socket, 
										  void (*cb)(void *),
										  void *prm);

/** puts socket in listening mode */
CPROPS_DLL
int cp_httpsocket_listen(cp_httpsocket *sock);

/** 
 * create a new cp_http_service descriptor. client requests for uris beginning
 * with 'path' will be delegated to this service once registered on a listening
 * socket.
 */
CPROPS_DLL
cp_http_service *cp_http_service_create(char *name, 
								        char *path, 
								        cp_http_service_callback service);
/** deallocate a cp_http_service descriptor */
CPROPS_DLL
void cp_http_service_delete(cp_http_service *svc);

/** register a service on a socket */
CPROPS_DLL
int cp_httpsocket_register_service(cp_httpsocket *server, cp_http_service *service);
/** unregister a service on a socket */
CPROPS_DLL
void *cp_httpsocket_unregister_service(cp_httpsocket *server, cp_http_service *service);

/** called when shutting down http layer */
CPROPS_DLL
void *cp_http_add_shutdown_callback(void (*cb)(void *), void *prm);

__END_DECLS
/** @} */

#endif

