//
// hooks_tpp.h
//
// instatiate at program start to handle install dir issues etc
//
// Copyright (c) 2006 Insilicos, LLC
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// $Author: bpratt $
//

#ifndef HOOKS_TPP_H
#define HOOKS_TPP_H

//
// instatiate at program start to handle install dir issues etc
//
class hooks_tpp { 
public:
	hooks_tpp(int &argc, char *argv[]); 
};

#endif // HOOKS_TPP_H
