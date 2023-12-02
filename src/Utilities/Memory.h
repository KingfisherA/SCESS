#pragma once
#include<pch.h>

template<typename A> using sRef = std::shared_ptr<A>;
template<typename A> using wRef = std::weak_ptr<A>;
template<typename A> using uRef = std::unique_ptr<A>;

template<typename A, typename... Types>
inline sRef<A> make(Types&&... arguments)
{
	return std::make_shared<A>(std::forward<Types>(arguments)...);
}

template<typename A, typename... Types>
inline uRef<A> umake(Types&&... arguments)
{
	return std::make_unique<A>(std::forward<Types>(arguments)...);
}