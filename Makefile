demo:
	g++ -Wall -O3 -Wextra -pedantic -o demo example.cpp ElpMpp02.cpp
trim_demo:
	g++ -Wall -O3 -Wextra -pedantic -o trim_demo example_usingElpMpp_trim.cpp ElpMpp02.cpp
javascript_demo:
	g++ -Wall -O3 -Wextra -pedantic -o javascript_demo example_usingElpMpp_JavaScript.cpp ElpMpp02.cpp
