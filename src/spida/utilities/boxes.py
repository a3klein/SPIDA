### Github CoPilot generated code to generate random boxes within a specified frame.
import random
import numpy as np
from typing import Tuple, Optional

def generate_random_box(frame_width: int, frame_height: int, 
                       box_width: int, box_height: int,
                       min_x: int = 0, min_y: int = 0,
                       max_x: Optional[int] = None, max_y: Optional[int] = None,
                       seed: Optional[int] = None) -> Tuple[int, int, int, int]:
    """
    Generate a random box within a specified frame.
    
    Parameters:
    -----------
    frame_width : int
        Width of the containing frame
    frame_height : int
        Height of the containing frame
    box_width : int
        Width of the box to generate
    box_height : int
        Height of the box to generate
    min_x : int, optional
        Minimum x coordinate for the box's top-left corner (default: 0)
    min_y : int, optional
        Minimum y coordinate for the box's top-left corner (default: 0)
    max_x : int, optional
        Maximum x coordinate for the box's top-left corner (default: frame_width - box_width)
    max_y : int, optional
        Maximum y coordinate for the box's top-left corner (default: frame_height - box_height)
    seed : int, optional
        Random seed for reproducibility
        
    Returns:
    --------
    Tuple[int, int, int, int]
        (x, y, width, height) where (x, y) is the top-left corner
        
    Raises:
    -------
    ValueError
        If the box dimensions are larger than the frame or invalid parameters
    """
    
    if seed is not None:
        random.seed(seed)
    
    # Validate inputs
    if box_width <= 0 or box_height <= 0:
        raise ValueError("Box dimensions must be positive")
    
    if frame_width <= 0 or frame_height <= 0:
        raise ValueError("Frame dimensions must be positive")
    
    if box_width > frame_width or box_height > frame_height:
        raise ValueError("Box dimensions cannot be larger than frame dimensions")
    
    # Set default max values if not provided
    if max_x is None:
        max_x = frame_width - box_width
    if max_y is None:
        max_y = frame_height - box_height
    
    # Validate coordinate constraints
    if min_x < 0 or min_y < 0:
        raise ValueError("Minimum coordinates must be non-negative")
    
    if max_x < min_x or max_y < min_y:
        raise ValueError("Maximum coordinates must be >= minimum coordinates")
    
    if min_x + box_width > frame_width or min_y + box_height > frame_height:
        raise ValueError("Box would exceed frame boundaries with given constraints")
    
    # Generate random coordinates
    x = random.randint(min_x, max_x)
    y = random.randint(min_y, max_y)
    
    return (x, y, box_width, box_height)


def generate_random_box_bounds(frame_width: int, frame_height: int, 
                              box_width: int, box_height: int,
                              **kwargs) -> Tuple[int, int, int, int]:
    """
    Generate a random box and return as (x1, y1, x2, y2) bounds.
    
    Returns:
    --------
    Tuple[int, int, int, int]
        (x1, y1, x2, y2) where (x1, y1) is top-left and (x2, y2) is bottom-right
    """
    x, y, w, h = generate_random_box(frame_width, frame_height, box_width, box_height, **kwargs)
    return (x, y, x + w, y + h)


def generate_random_box_center(frame_width: int, frame_height: int, 
                              box_width: int, box_height: int,
                              **kwargs) -> Tuple[int, int, int, int]:
    """
    Generate a random box and return with center coordinates.
    
    Returns:
    --------
    Tuple[int, int, int, int]
        (center_x, center_y, width, height)
    """
    x, y, w, h = generate_random_box(frame_width, frame_height, box_width, box_height, **kwargs)
    center_x = x + w // 2
    center_y = y + h // 2
    return (center_x, center_y, w, h)


def generate_multiple_random_boxes(frame_width: int, frame_height: int,
                                  box_width: int, box_height: int,
                                  num_boxes: int,
                                  avoid_overlap: bool = False,
                                  max_attempts: int = 1000,
                                  **kwargs) -> list:
    """
    Generate multiple random boxes within a frame.
    
    Parameters:
    -----------
    num_boxes : int
        Number of boxes to generate
    avoid_overlap : bool
        If True, ensure boxes don't overlap (may fail if too many boxes requested)
    max_attempts : int
        Maximum attempts to place non-overlapping boxes
        
    Returns:
    --------
    list
        List of (x, y, width, height) tuples
    """
    boxes = []
    
    if not avoid_overlap:
        # Simple case: just generate random boxes
        for _ in range(num_boxes):
            box = generate_random_box(frame_width, frame_height, box_width, box_height, **kwargs)
            boxes.append(box)
    else:
        # Complex case: avoid overlaps
        attempts = 0
        while len(boxes) < num_boxes and attempts < max_attempts:
            candidate = generate_random_box(frame_width, frame_height, box_width, box_height, **kwargs)
            
            # Check for overlap with existing boxes
            overlaps = False
            for existing_box in boxes:
                if boxes_overlap(candidate, existing_box):
                    overlaps = True
                    break
            
            if not overlaps:
                boxes.append(candidate)
            
            attempts += 1
        
        if len(boxes) < num_boxes:
            print(f"Warning: Could only place {len(boxes)} out of {num_boxes} non-overlapping boxes")
    
    return boxes


def boxes_overlap(box1: Tuple[int, int, int, int], box2: Tuple[int, int, int, int]) -> bool:
    """
    Check if two boxes overlap.
    
    Parameters:
    -----------
    box1, box2 : Tuple[int, int, int, int]
        Boxes in (x, y, width, height) format
        
    Returns:
    --------
    bool
        True if boxes overlap, False otherwise
    """
    x1, y1, w1, h1 = box1
    x2, y2, w2, h2 = box2
    
    # Convert to bounds
    x1_end, y1_end = x1 + w1, y1 + h1
    x2_end, y2_end = x2 + w2, y2 + h2
    
    # Check for no overlap conditions
    if x1_end <= x2 or x2_end <= x1 or y1_end <= y2 or y2_end <= y1:
        return False
    
    return True


def visualize_boxes(frame_width: int, frame_height: int, boxes: list, ax=None, 
                   show_plot: bool = True, save_path: Optional[str] = None):
    """
    Visualize the generated boxes within the frame.
    
    Parameters:
    -----------
    boxes : list
        List of (x, y, width, height) tuples
    show_plot : bool
        Whether to display the plot
    save_path : str, optional
        Path to save the plot
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        
        if ax is None: 
            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        # Draw frame
        frame_rect = patches.Rectangle((0, 0), frame_width, frame_height, 
                                     linewidth=2, edgecolor='black', facecolor='lightgray', alpha=0.3)
        ax.add_patch(frame_rect)
        
        # Draw boxes
        colors = plt.cm.Set3(np.linspace(0, 1, len(boxes)))
        for i, (x, y, w, h) in enumerate(boxes):
            rect = patches.Rectangle((x, y), w, h, 
                                   linewidth=2, edgecolor='red', facecolor=colors[i], alpha=0.7)
            ax.add_patch(rect)
            
            # Add box number
            ax.text(x + w/2, y + h/2, str(i+1), ha='center', va='center', 
                   fontsize=12, fontweight='bold')
        
        ax.set_xlim(-10, frame_width + 10)
        ax.set_ylim(-10, frame_height + 10)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Random Boxes in {frame_width}x{frame_height} Frame')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        if show_plot:
            plt.show()
            
    except ImportError:
        print("Matplotlib not available for visualization")