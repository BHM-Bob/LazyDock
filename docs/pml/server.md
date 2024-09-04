# Overview
This module provides a client-server architecture for interacting with PyMOL, a molecular visualization system. The `VServer` class runs on the server side, listening for connections and executing PyMOL commands sent by the client. The `VClient` class is used on the client side to send commands to the server. The `PymolAPI` class acts as a wrapper for PyMOL's API, providing a convenient interface for command execution.

# Classes
## `PymolAPI`
### Description
A wrapper class for PyMOL's API, allowing for the execution of PyMOL commands (`cmd`) or API functions (`api`) through a client-server architecture.

### Methods
#### `send_action(api: str, fn: str, *args, **kwargs) -> Any`
##### Description
Sends an action to the server to be executed.

##### Parameters
- `api`: The API type ('cmd' or 'api').
- `fn`: The function name to execute.
- `*args`: Positional arguments for the function.
- `**kwargs`: Keyword arguments for the function.

##### Returns
- `Any`: The return value of the executed function.

## `VServer`
### Description
A server class that listens for connections from clients and executes PyMOL commands sent by the clients.

### Methods
#### `copy_logs() -> List[str]`
##### Description
Copies the server logs.

##### Returns
- `List[str]`: A list of log strings.

#### `copy_new_logs() -> List[str]`
##### Description
Copies new logs since the last check.

##### Returns
- `List[str]`: A list of new log strings.

#### `_add_log(log: str, verbose: bool = False) -> None`
##### Description
Adds a log entry.

##### Parameters
- `log`: The log message.
- `verbose`: Whether to print the log message.

#### `_recvall(sock: socket.socket, count: int)`
##### Description
Receives all data from the socket.

##### Parameters
- `sock`: The socket object.
- `count`: The number of bytes to receive.

#### `_sendall(scok: socket.socket, data: bytes) -> None`
##### Description
Sends all data to the socket.

##### Parameters
- `scok`: The socket object.
- `data`: The data to send.

#### `_run_server() -> None`
##### Description
The main loop of the server, running in a separate thread, listening for connections and executing commands.

## `VClient`
### Description
A client class that connects to the server and sends commands to be executed in PyMOL.

### Methods
#### `send_action(api: str, fn: str, *args, **kwargs) -> Any`
##### Description
Sends an action to the server to be executed.

##### Parameters
- `api`: The API type ('cmd' or 'api').
- `fn`: The function name to execute.
- `*args`: Positional arguments for the function.
- `**kwargs`: Keyword arguments for the function.

##### Returns
- `Any`: The return value of the executed function.

#### `sen_quit() -> None`
##### Description
Sends a quit command to the server and closes the connection.

## Functions
### `_test_server()`
##### Description
A test function to start the server and print logs.

### `_test_client()`
##### Description
A test function to start the client, connect to the server, and execute a PyMOL command.

# Main Execution
The provided code includes a main execution block that starts both the server and client in separate processes, allowing for testing of the client-server architecture.
