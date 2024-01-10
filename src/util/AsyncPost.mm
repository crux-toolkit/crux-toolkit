// AsyncPost.mm
#import <Foundation/Foundation.h>

extern "C" void performAsyncPOSTRequest(const char *urlString, const char *jsonString) {
    @autoreleasepool {
        // Convert C strings to NSStrings
        NSString *urlStringObjC = [NSString stringWithUTF8String:urlString];
        NSString *jsonStringObjC = [NSString stringWithUTF8String:jsonString];

        // Create a URL
        NSURL *url = [NSURL URLWithString:urlStringObjC];

        // Create a mutable request
        NSMutableURLRequest *request = [NSMutableURLRequest requestWithURL:url];
        [request setHTTPMethod:@"POST"];

        // Set Content-Type header
        [request setValue:@"application/json" forHTTPHeaderField:@"Content-Type"];

        // Set User-Agent header
        [request setValue:@"Crux" forHTTPHeaderField:@"User-Agent"];

        // Convert JSON string to NSData
        NSData *jsonData = [jsonStringObjC dataUsingEncoding:NSUTF8StringEncoding];
        [request setHTTPBody:jsonData];

        // Create a session configuration
        NSURLSessionConfiguration *config = 
          [NSURLSessionConfiguration defaultSessionConfiguration];

        // Create a session
        NSURLSession *session = [NSURLSession sessionWithConfiguration:config];

        // Create a task with a completion handler
        NSURLSessionDataTask *task = 
          [session dataTaskWithRequest:request completionHandler:^(NSData * _Nullable data, NSURLResponse * _Nullable response, NSError * _Nullable error) {
            if (error) {
                fprintf(stderr, "Error: %s\n", [error.localizedDescription UTF8String]);
            } else {
              // Get the HTTP response code
              NSHTTPURLResponse *httpResponse = (NSHTTPURLResponse *)response;
              NSInteger statusCode = httpResponse.statusCode;
              // Print the HTTP response code
              // printf("HTTP Response Code: %ld\n", (long)statusCode);

            }
            // Stop the run loop to exit the program
            CFRunLoopStop(CFRunLoopGetCurrent());
        }];

        // Perform the task
        [task resume];

        // Run the run loop to keep the program alive while waiting for the asynchronous task to complete
        // CFRunLoopRun();
    }
}

