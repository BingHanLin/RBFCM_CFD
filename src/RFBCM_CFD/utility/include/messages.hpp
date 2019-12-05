#ifndef MESSAGES_HPP
#define MESSAGES_HPP
// https://stackoverflow.com/questions/2193544/how-to-print-additional-information-when-assert-fails
#define ASSERT(condition)                                                \
    {                                                                    \
        std::cerr << "ASSERT FAILED: " << condition << " @ " << __FILE__ \
                  << " ( line: " << __LINE__ << ")" << std::endl;        \
        exit(-1);                                                        \
    }

#endif
